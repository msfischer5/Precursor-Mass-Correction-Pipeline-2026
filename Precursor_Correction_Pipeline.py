"""
UniDec (MS1) + TopFD msalign (MS2) -> Corrected msalign for TopPIC

Steps
1. Load mzML; record the precursor window and parent MS1 scan of every MS2
2. Extract the parent MS1 scans (only those that precede MS2 activity)
3. Deconvolute each MS1 scan with UniDec and build a mass list
4. Parse the TopFD ms2.msalign file
5. For each MS2, search the UniDec masses of its parent MS1 scan for a match
   to the precursor window midpoint within the given PPM tolerance
6. Apply match constraints (mass threshold, noisy-scan rejection,
   RT-replication corroboration)
7. Write a new msalign file containing the new precursor masses and the
   original fragments (replaced spectra only)

Run `python unidec_precursor_correction.py --help` for command-line usage.
"""

import argparse
import sys
import io
import os
import re
import shutil
import numpy as np
import pyopenms

from collections import Counter
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from unidec.engine import UniDec

PROTON_MASS = 1.007276

# Name of the file written into each output directory.
OUTPUT_FILENAME = "unidec_corrected_ms2.msalign"

DEFAULTS = dict(
    mass_threshold=30000.0,
    mass_ppm_tolerance=100.0,
    match_charge_range=(5, 115),
    unidec_mass_range=(4000, 110000),
    unidec_charge_range=(5, 115),
    peak_thresh=0.1,
    peak_window=15.0,
    max_peaks_per_scan=15,
    n_workers=None,
    rt_replication_window_min=1.0,
    min_rt_occurrences=2,
    background_subtype=2,     # 0 = off, 1 = line, 2 = curved (may differ by UniDec version)
    background_subbuff=100,
)

# UniDec settings that are fixed unless overridden.
DEFAULT_UNIDEC_SETTINGS = dict(
    minmz=500,        # m/z axis lower bound
    maxmz=3000,       # m/z axis upper bound
    mzbins=1.0,       # m/z bin size
    massbins=1.0,     # mass bin size (Da)
    psig=10.0,        # peak shape FWHM in mass-space smoothing
    psfun=0,          # peak shape function (0 = Gaussian)
    beta=100,         # charge-state smoothing / artifact suppression
    mzsig=0.5,        # m/z peak FWHM
)


def _deconvolve_one_scan_worker(scan_num, txt_path_str, mass_lb, mass_ub,
                                 start_z, end_z, peak_thresh, peak_window,
                                 max_peaks_per_scan,
                                 background_subtype=0, background_subbuff=0,
                                 unidec_settings=None):
    """
    Run UniDec on a single MS1 scan with multiprocessing support.

    Args:
        scan_num:            MS1 scan number
        txt_path_str:        Path to the two-column m/z/intensity txt
        mass_lb, mass_ub:    Mass axis bounds for this scan's UniDec run (Da)
        start_z, end_z:      Charge range UniDec considers
        peak_thresh:         UniDec peak-picker intensity threshold
        peak_window:         UniDec peak-picker mass separation (Da)
        max_peaks_per_scan:  Maximum number of picked peaks before scan rejection
        background_subtype:  Controls the type of UniDec background subtraction used
        background_subbuff:  Background subtraction parameter
        unidec_settings:     Dict of additional UniDec config attributes,
                             merged over DEFAULT_UNIDEC_SETTINGS

    Returns:
        (scan_num, masses, was_noisy, error_message_or_None)
        masses is a list of (mass, intensity) tuples, empty if rejected/failed.
    """
    txt_path = Path(txt_path_str)
    old_out, old_err = sys.stdout, sys.stderr
    settings = {**DEFAULT_UNIDEC_SETTINGS, **(unidec_settings or {})}
    try:
        sys.stdout = sys.stderr = io.StringIO()

        ud = UniDec()
        ud.open_file(str(txt_path))

        for key, value in settings.items():
            setattr(ud.config, key, value)

        ud.config.masslb     = mass_lb
        ud.config.massub     = mass_ub
        ud.config.startz     = start_z
        ud.config.endz       = end_z
        ud.config.peakthresh = peak_thresh
        ud.config.peakwindow = peak_window
        ud.config.subtype    = background_subtype
        ud.config.subbuff    = background_subbuff

        ud.process_data()
        ud.run_unidec()
        ud.pick_peaks()

        sys.stdout, sys.stderr = old_out, old_err

        masses = []
        pks = getattr(ud, 'pks', None)
        if pks is not None and getattr(pks, 'peaks', None):
            for p in pks.peaks:
                mass = getattr(p, 'mass', None)
                if mass is None:
                    continue
                intensity = getattr(p, 'height', None)
                if intensity is None:
                    intensity = getattr(p, 'integral', None)
                if intensity is None:
                    continue
                masses.append((float(mass), float(intensity)))

        if len(masses) > max_peaks_per_scan:
            return (scan_num, [], True, None)
        return (scan_num, masses, False, None)

    except Exception as e:
        sys.stdout, sys.stderr = old_out, old_err
        return (scan_num, [], False, str(e))

    finally:
        unidec_folder = txt_path.parent / f"{txt_path.stem}_unidecfiles"
        shutil.rmtree(unidec_folder, ignore_errors=True)
        txt_path.unlink(missing_ok=True)


class UniDecDDAProcessor:

    def __init__(self, mzml_path, output_dir, topfd_msalign_path):
        """
        Args:
            mzml_path:          Path to input mzML file
            output_dir:         Directory for all outputs
            topfd_msalign_path: Path to TopFD ms2.msalign file
        """
        self.mzml_path          = Path(mzml_path)
        self.output_dir         = Path(output_dir)
        self.topfd_msalign_path = Path(topfd_msalign_path)

        self.ms1_dir = self.output_dir / "MS1_scans"
        self.ms1_dir.mkdir(parents=True, exist_ok=True)

        self.ms1_scan_map: dict = {}
        self.ms2_precursor_map: dict = {}
        self.unidec_results: dict = {}
        self.noisy_scans: set = set()


    def load_mzml(self):
        print(f"Loading mzML: {self.mzml_path}")
        self.exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(str(self.mzml_path), self.exp)
        total = self.exp.size()

        level_counts = Counter(int(s.getMSLevel()) for s in self.exp)
        level_str = "  ".join(f"MS{k}: {v}" for k, v in sorted(level_counts.items()))
        print(f"  Loaded {total} spectra  ({level_str})")

    @staticmethod
    def _native_scan_number(spectrum, fallback: int) -> int:
        """Parse the integer scan number from a spectrum's native ID string."""
        native_id = spectrum.getNativeID()
        for token in ('scan', 'index'):
            m = re.search(rf'{token}=(\d+)', native_id)
            if m:
                return int(m.group(1))
        return fallback


    def extract_ms1_scans(self, needed_ms1_scans: set):
        """
        Save MS1 scans as two-column txt files for UniDec.

        Only scans whose native scan number appears in needed_ms1_scans are
        written.  This set is built from extract_ms2_precursors() so that
        UniDec only runs on MS1 scans that precede at least one MS2
        spectrum.

        Args:
            needed_ms1_scans: Set of native MS1 scan numbers to process.
        """
        print(f"\nExtracting MS1 scans (only those preceding MS2 activity)...")
        count = 0
        skipped = 0
        for i, spectrum in enumerate(self.exp):
            if int(spectrum.getMSLevel()) != 1:
                continue
            scan_num = self._native_scan_number(spectrum, i + 1)
            if scan_num not in needed_ms1_scans:
                skipped += 1
                continue
            mz, intensity = spectrum.get_peaks()
            if len(mz) < 10:
                continue
            txt_path = self.ms1_dir / f"scan_{scan_num}.txt"
            np.savetxt(txt_path, np.column_stack((mz, intensity)),
                       fmt='%.6f %.6f', delimiter=' ')
            self.ms1_scan_map[scan_num] = txt_path
            count += 1
        print(f"  Extracted {count} MS1 scans (skipped {skipped} with no MS2 activity)")


    def extract_ms2_precursors(self):
        """
        Record the precursor window midpoint m/z and parent MS1 scan for
        every MS2 spectrum.

        """
        print("\nExtracting MS2 precursor windows...")
        count = 0
        for i, spectrum in enumerate(self.exp):
            if int(spectrum.getMSLevel()) != 2:
                continue
            precursors = spectrum.getPrecursors()
            if not precursors:
                continue

            p         = precursors[0]
            scan_num  = self._native_scan_number(spectrum, i + 1)
            center_mz = p.getMZ()
            win_lower = p.getIsolationWindowLowerOffset()
            win_upper = p.getIsolationWindowUpperOffset()

            if win_lower == 0.0 and win_upper == 0.0:
                midpoint = center_mz
            else:
                midpoint = ((center_mz - win_lower) + (center_mz + win_upper)) / 2.0

            self.ms2_precursor_map[scan_num] = {
                'midpoint_mz': midpoint,
                'ms1_scan':    self._find_preceding_ms1(i),
                'rt':          spectrum.getRT(),
            }
            count += 1
        print(f"  Extracted window info for {count} MS2 spectra")

    def _find_preceding_ms1(self, ms2_index: int):
        """Return native scan number of the MS1 scan just before ms2_index."""
        for i in range(ms2_index - 1, -1, -1):
            if int(self.exp[i].getMSLevel()) == 1:
                return self._native_scan_number(self.exp[i], i + 1)
        return None


    def run_unidec_batch(self,
                          mass_range=DEFAULTS['unidec_mass_range'],
                          charge_range=DEFAULTS['unidec_charge_range'],
                          peak_thresh=DEFAULTS['peak_thresh'],
                          peak_window=DEFAULTS['peak_window'],
                          max_peaks_per_scan=DEFAULTS['max_peaks_per_scan'],
                          n_workers=DEFAULTS['n_workers'],
                          background_subtype=DEFAULTS['background_subtype'],
                          background_subbuff=DEFAULTS['background_subbuff'],
                          unidec_settings=None):
        """
        Run UniDec on every extracted MS1 scan in parallel.

        Args:
            mass_range:          (min_mass, max_mass) in Da
            charge_range:        (min_z, max_z) UniDec itself searches
            peak_thresh:         UniDec config.peakthresh — minimum peak
                                  height as a fraction of the tallest peak in
                                  the mass spectrum (0-1).
            peak_window:         UniDec config.peakwindow (Da) — minimum mass
                                  separation between two calls being counted
                                  as distinct peaks.
            max_peaks_per_scan:  If UniDec's picker finds more peaks than
                                  this in a scan, the scan is treated as too
                                  noisy/unresolved to trust and all of its
                                  masses are discarded.
            n_workers:           Number of parallel worker processes
                                  (default: CPU count - 1).
            background_subtype:  UniDec config.subtype — background
                                  subtraction mode.
            background_subbuff:  UniDec config.subbuff — parameter that
                                  accompanies background_subtype.
            unidec_settings:     Optional dict overriding entries of
                                  DEFAULT_UNIDEC_SETTINGS.
        """
        print("\nRunning UniDec on MS1 scans...")
        total = len(self.ms1_scan_map)

        if n_workers is None:
            n_workers = max(1, (os.cpu_count() or 2) - 1)
        print(f"  Using {n_workers} parallel worker process(es) "
              f"(detected {os.cpu_count()} CPU(s))")

        jobs = [
            (scan_num, str(txt_path), mass_range[0], mass_range[1],
             charge_range[0], charge_range[1],
             peak_thresh, peak_window, max_peaks_per_scan,
             background_subtype, background_subbuff, unidec_settings)
            for scan_num, txt_path in sorted(self.ms1_scan_map.items())
        ]

        noisy_count = 0
        completed = 0
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(_deconvolve_one_scan_worker, *job): job[0]
                       for job in jobs}
            for future in as_completed(futures):
                scan_num = futures[future]
                try:
                    s, masses, was_noisy, err = future.result()
                except Exception as e:
                    s, masses, was_noisy, err = scan_num, [], False, str(e)

                completed += 1
                if err:
                    print(f"  Scan {s} ({completed}/{total}): ERROR: {err}")
                    self.unidec_results[s] = []
                elif was_noisy:
                    self.noisy_scans.add(s)
                    noisy_count += 1
                    self.unidec_results[s] = []
                    print(f"  Scan {s} ({completed}/{total}): "
                          f"REJECTED (noisy, > {max_peaks_per_scan} peaks)")
                else:
                    self.unidec_results[s] = masses
                    print(f"  Scan {s} ({completed}/{total}): {len(masses)} peaks")

        n_ok = sum(1 for v in self.unidec_results.values() if v)
        print(f"\n  {n_ok}/{len(self.unidec_results)} MS1 scans yielded usable UniDec peaks")
        print(f"  {noisy_count} scan(s) rejected as noisy (> {max_peaks_per_scan} peaks)")


    @staticmethod
    def _is_fragment_line(stripped: str) -> bool:
        """
        """
        if not stripped or '=' in stripped:
            return False
        parts = stripped.split()
        if len(parts) < 2:
            return False
        try:
            float(parts[0])
            float(parts[1])
            return True
        except ValueError:
            return False

    def parse_topfd_msalign(self):
        """
        Parse the TopFD ms2.msalign file.

        Returns:
            dict: native scan_num (int) -> entry dict with keys:
                  'header_lines'        – list of raw str (all non-fragment lines)
                  'fragment_lines'      – list of raw str
                  'scan_number'         – int, same as key
                  'expected_fragments'  – int from DECONVOLUTED_MASS_NUMBER or None
                  'win_begin'           – float or None
                  'win_end'             – float or None
        """
        print(f"\nParsing TopFD msalign: {self.topfd_msalign_path}")

        if not self.topfd_msalign_path.exists():
            raise FileNotFoundError(
                f"TopFD msalign not found: {self.topfd_msalign_path}"
            )

        entries            = {}
        current_entry      = None
        current_scan       = None
        fragment_lines     = []
        mismatch_warnings  = 0

        with open(self.topfd_msalign_path, 'r') as fh:
            for raw_line in fh:
                stripped = raw_line.strip()

                if stripped == 'BEGIN IONS':
                    current_entry = {
                        'header_lines':      [],
                        'fragment_lines':    [],
                        'scan_number':       None,
                        'expected_fragments': None,
                        'win_begin':         None,
                        'win_end':           None,
                    }
                    fragment_lines = []

                elif stripped.startswith('SCANS='):
                    current_scan = int(stripped.split('=')[1])
                    current_entry['scan_number'] = current_scan
                    current_entry['header_lines'].append(raw_line)

                elif stripped.startswith('PRECURSOR_WINDOW_BEGIN='):
                    current_entry['win_begin'] = float(stripped.split('=')[1])
                    current_entry['header_lines'].append(raw_line)

                elif stripped.startswith('PRECURSOR_WINDOW_END='):
                    current_entry['win_end'] = float(stripped.split('=')[1])
                    current_entry['header_lines'].append(raw_line)

                elif stripped.startswith('DECONVOLUTED_MASS_NUMBER='):
                    try:
                        current_entry['expected_fragments'] = int(
                            stripped.split('=')[1].strip()
                        )
                    except ValueError:
                        pass
                    current_entry['header_lines'].append(raw_line)

                elif stripped == 'END IONS':
                    if current_entry is not None and current_scan is not None:
                        current_entry['fragment_lines'] = fragment_lines
                        expected = current_entry['expected_fragments']
                        captured = len(fragment_lines)
                        if expected is not None and captured != expected:
                            print(f"  WARNING: scan {current_scan} – "
                                  f"DECONVOLUTED_MASS_NUMBER={expected} "
                                  f"but captured {captured} fragment lines")
                            mismatch_warnings += 1
                        entries[current_scan] = current_entry
                    current_entry  = None
                    current_scan   = None
                    fragment_lines = []

                elif current_entry is not None:
                    if self._is_fragment_line(stripped):
                        fragment_lines.append(raw_line)
                    else:
                        current_entry['header_lines'].append(raw_line)

        print(f"  Parsed {len(entries)} msalign entries")
        if mismatch_warnings:
            print(f"  WARNING: {mismatch_warnings} entries had fragment count mismatches")
        return entries


    def _find_best_mass_match(self, query_mz: float, unidec_masses: list,
                               mass_ppm_tolerance: float, charge_range: tuple):
        """
        Search every charge state in charge_range. For each z, compute the
        implied neutral mass and see if it lands within mass_ppm_tolerance
        of any UniDec-picked mass. Returns the highest-intensity match
        found across all tested charges, or None if nothing matched.

        Returns:
            (mass, charge, intensity) or None
        """
        best_match = None
        best_intensity = -1.0

        for z in range(charge_range[0], charge_range[1] + 1):
            query_mass = query_mz * z - z * PROTON_MASS
            for mass, intensity in unidec_masses:
                mass_ppm_error = abs(mass - query_mass) / query_mass * 1e6
                if mass_ppm_error < mass_ppm_tolerance:
                    if intensity > best_intensity:
                        best_intensity = intensity
                        best_match = (mass, z, intensity)

        return best_match

    @staticmethod
    def _apply_replication_filter(candidates: dict, mass_ppm_tolerance: float,
                                   rt_replication_window_min: float):
        """
        Given raw candidate replacements (scan_num -> dict with 'mass' and
        'rt' in seconds), count, for each candidate, how many other
        candidates have a closely-matching mass (within mass_ppm_tolerance)
        within rt_replication_window_min minutes.
        """
        scan_nums = list(candidates.keys())
        window_sec = rt_replication_window_min * 60.0

        for scan_i in scan_nums:
            cand_i = candidates[scan_i]
            occurrences = 1  # counts itself
            for scan_j in scan_nums:
                if scan_j == scan_i:
                    continue
                cand_j = candidates[scan_j]
                mass_ppm_error = (abs(cand_i['mass'] - cand_j['mass'])
                                  / cand_j['mass'] * 1e6)
                if mass_ppm_error > mass_ppm_tolerance:
                    continue
                if abs(cand_i['rt'] - cand_j['rt']) > window_sec:
                    continue
                occurrences += 1
            cand_i['rt_occurrences'] = occurrences

        return candidates

    def find_unidec_replacements(self, msalign_entries: dict,
                                  mass_ppm_tolerance: float = DEFAULTS['mass_ppm_tolerance'],
                                  charge_range: tuple = DEFAULTS['match_charge_range'],
                                  mass_threshold: float = DEFAULTS['mass_threshold'],
                                  rt_replication_window_min: float = DEFAULTS['rt_replication_window_min'],
                                  min_rt_occurrences: int = DEFAULTS['min_rt_occurrences']) -> dict:
        """
        For each MS2, compute the precursor window midpoint and attempt to
        match it against UniDec masses from the preceding MS1 scan.

        Replacement requires:
          a) A UniDec mass matches the midpoint within mass_ppm_tolerance
          b) That mass >= mass_threshold
          c) Corroboration: the matched mass is seen again (within
             mass_ppm_tolerance) for another MS2 spectrum within
             rt_replication_window_min, for a total of >= min_rt_occurrences
             occurrences (including itself).

        Args:
            msalign_entries:            Output of parse_topfd_msalign()
            mass_ppm_tolerance:         PPM tolerance for midpoint matching
                                         (also used for replication matching)
            charge_range:               (min_z, max_z) charge states to test
            mass_threshold:             Minimum UniDec mass (Da) for replacement
            rt_replication_window_min:  RT window (minutes) for replication check
            min_rt_occurrences:         Occurrences (including itself) required

        Returns:
            dict: scan_num -> {mass, charge, mz, intensity, ms1_scan, rt}
                  Contains only spectra where replacement was made.
        """
        print(f"\nSearching for UniDec replacements...")
        print(f"  Mass PPM tolerance: {mass_ppm_tolerance}")
        print(f"  Mass threshold: {mass_threshold/1000:.0f} kDa")
        print(f"  RT-replication window: {rt_replication_window_min} min "
              f"(need {min_rt_occurrences} occurrences)")

        candidates = {}
        stats = {'replaced': 0, 'below_threshold': 0,
                 'no_match': 0, 'no_ms1': 0, 'uncorroborated': 0}

        for scan_num, entry in msalign_entries.items():

            win_begin = entry.get('win_begin')
            win_end   = entry.get('win_end')

            if win_begin is not None and win_end is not None:
                midpoint = (win_begin + win_end) / 2.0
            elif scan_num in self.ms2_precursor_map:
                # Fall back to mzML-derived midpoint
                midpoint = self.ms2_precursor_map[scan_num]['midpoint_mz']
            else:
                stats['no_ms1'] += 1
                continue

            ms1_scan = None
            rt       = 0.0
            if scan_num in self.ms2_precursor_map:
                ms1_scan = self.ms2_precursor_map[scan_num]['ms1_scan']
                rt       = self.ms2_precursor_map[scan_num]['rt']

            if ms1_scan is None or ms1_scan not in self.unidec_results:
                stats['no_ms1'] += 1
                continue

            unidec_masses = self.unidec_results[ms1_scan]
            if not unidec_masses:
                stats['no_ms1'] += 1
                continue

            match = self._find_best_mass_match(
                midpoint, unidec_masses, mass_ppm_tolerance, charge_range,
            )

            if match is None:
                stats['no_match'] += 1
                continue

            unidec_mass, unidec_charge, unidec_intensity = match

            if unidec_mass < mass_threshold:
                stats['below_threshold'] += 1
                continue

            candidates[scan_num] = {
                'mass':      unidec_mass,
                'charge':    unidec_charge,
                'mz':        midpoint,
                'intensity': unidec_intensity,
                'ms1_scan':  ms1_scan,
                'rt':        rt,
            }

        candidates = self._apply_replication_filter(
            candidates, mass_ppm_tolerance, rt_replication_window_min,
        )

        replacements = {}
        for scan_num, cand in candidates.items():
            if cand['rt_occurrences'] < min_rt_occurrences:
                stats['uncorroborated'] += 1
                continue

            replacements[scan_num] = {
                'mass':      cand['mass'],
                'charge':    cand['charge'],
                'mz':        cand['mz'],
                'intensity': cand['intensity'],
                'ms1_scan':  cand['ms1_scan'],
                'rt':        cand['rt'],
            }
            stats['replaced'] += 1

        print(f"  Replaced           : {stats['replaced']}")
        print(f"  Below threshold    : {stats['below_threshold']}")
        print(f"  No UniDec match    : {stats['no_match']}")
        print(f"  No MS1 / UniDec    : {stats['no_ms1']}")
        print(f"  Uncorroborated     : {stats['uncorroborated']} "
              f"(matched but not re-isolated nearby in RT)")
        return replacements

    # Output
    def write_corrected_msalign(self, replacements: dict,
                                 msalign_entries: dict,
                                 output_file=None):
        """
        Write a new msalign containing the replaced spectra.

        For each replaced spectrum every header line is copied verbatim
        except the five precursor fields below, which are rewritten:

            PRECURSOR_MZ        = window midpoint
            PRECURSOR_CHARGE    = UniDec matched charge
            PRECURSOR_MASS      = UniDec matched mass
            PRECURSOR_INTENSITY = UniDec peak intensity
            PRECURSOR_FEATURE_ID = -1

        Args:
            replacements:    Output of find_unidec_replacements()
            msalign_entries: Output of parse_topfd_msalign()
            output_file:     Destination path
                             (default: output_dir/OUTPUT_FILENAME)
        """
        if output_file is None:
            output_file = self.output_dir / OUTPUT_FILENAME

        print(f"\nWriting corrected msalign -> {output_file}")
        print(f"  Spectra to write: {len(replacements)}")

        REPLACED_FIELDS = {
            'PRECURSOR_MZ=',
            'PRECURSOR_CHARGE=',
            'PRECURSOR_MASS=',
            'PRECURSOR_INTENSITY=',
            'PRECURSOR_FEATURE_ID=',
        }

        written = 0
        with open(output_file, 'w') as fh:
            fh.write('#Created by UniDec-TopFD Pipeline\n')
            fh.write(f'#Timestamp: {datetime.now().strftime("%a %b %d %H:%M:%S %Y")}\n')
            fh.write(f'#Source mzML: {self.mzml_path}\n')
            fh.write(f'#TopFD msalign: {self.topfd_msalign_path}\n')
            fh.write('#UniDec-corrected spectra only\n')
            fh.write(f'#Number of MS/MS spectra: {len(replacements)}\n')
            fh.write('\n')

            for scan_num in sorted(replacements.keys()):
                rep   = replacements[scan_num]
                entry = msalign_entries[scan_num]

                fh.write('BEGIN IONS\n')

                for line in entry['header_lines']:
                    # Identify the field key robustly (handles multi-value lines)
                    key = (line.split('=')[0] + '=') if '=' in line else ''
                    if key in REPLACED_FIELDS:
                        if key == 'PRECURSOR_MZ=':
                            fh.write(f'PRECURSOR_MZ={rep["mz"]:.6f}\n')
                        elif key == 'PRECURSOR_CHARGE=':
                            fh.write(f'PRECURSOR_CHARGE={rep["charge"]}\n')
                        elif key == 'PRECURSOR_MASS=':
                            fh.write(f'PRECURSOR_MASS={rep["mass"]:.5f}\n')
                        elif key == 'PRECURSOR_INTENSITY=':
                            fh.write(f'PRECURSOR_INTENSITY={rep["intensity"]:.2f}\n')
                        elif key == 'PRECURSOR_FEATURE_ID=':
                            fh.write('PRECURSOR_FEATURE_ID=-1\n')
                    else:
                        fh.write(line)

                for line in entry['fragment_lines']:
                    fh.write(line)

                fh.write('END IONS\n\n')
                written += 1

        print(f"  Wrote {written} spectra")
        return output_file


    def run_full_pipeline(self,
                          mass_threshold: float = DEFAULTS['mass_threshold'],
                          mass_ppm_tolerance: float = DEFAULTS['mass_ppm_tolerance'],
                          match_charge_range: tuple = DEFAULTS['match_charge_range'],
                          unidec_mass_range: tuple = DEFAULTS['unidec_mass_range'],
                          unidec_charge_range: tuple = DEFAULTS['unidec_charge_range'],
                          peak_thresh: float = DEFAULTS['peak_thresh'],
                          peak_window: float = DEFAULTS['peak_window'],
                          max_peaks_per_scan: int = DEFAULTS['max_peaks_per_scan'],
                          n_workers: int = DEFAULTS['n_workers'],
                          rt_replication_window_min: float = DEFAULTS['rt_replication_window_min'],
                          min_rt_occurrences: int = DEFAULTS['min_rt_occurrences'],
                          background_subtype: int = DEFAULTS['background_subtype'],
                          background_subbuff: float = DEFAULTS['background_subbuff'],
                          unidec_settings: dict = None):
        """
        Run the complete pipeline.

        Args:
            mass_threshold:             Minimum UniDec mass (Da) to trigger replacement
            mass_ppm_tolerance:         PPM tolerance for precursor midpoint matching
            match_charge_range:         Charge states tested during UniDec m/z search
            unidec_mass_range:          (min, max) mass window passed to UniDec (Da)
            unidec_charge_range:        (min, max) charge window passed to UniDec
            peak_thresh:                UniDec native peak-picker intensity threshold
            peak_window:                UniDec native peak-picker mass separation (Da)
            max_peaks_per_scan:         Reject any MS1 scan with more picked peaks than this.
            n_workers:                  Parallel worker processes for UniDec deconvolution.
                                         Defaults to (CPU count - 1).
            rt_replication_window_min:  RT window (min) for replication corroboration
            min_rt_occurrences:         Occurrences (including itself) required.
            background_subtype:         UniDec config.subtype — background subtraction
                                         mode.
            background_subbuff:         UniDec config.subbuff — parameter that
                                         accompanies background_subtype.
            unidec_settings:            Optional dict overriding entries of
                                         DEFAULT_UNIDEC_SETTINGS (e.g. {'beta': 50}).

        Returns:
            Path of the written msalign file.
        """
        print("=" * 60)
        print("DDA Top-Down Proteomics Pipeline")
        print("=" * 60)

        self.load_mzml()

        self.extract_ms2_precursors()
        needed_ms1_scans = {
            info['ms1_scan']
            for info in self.ms2_precursor_map.values()
            if info['ms1_scan'] is not None
        }
        print(f"  {len(needed_ms1_scans)} unique MS1 scans precede MS2 activity")

        self.extract_ms1_scans(needed_ms1_scans)
        self.run_unidec_batch(mass_range=unidec_mass_range,
                              charge_range=unidec_charge_range,
                              peak_thresh=peak_thresh,
                              peak_window=peak_window,
                              max_peaks_per_scan=max_peaks_per_scan,
                              n_workers=n_workers,
                              background_subtype=background_subtype,
                              background_subbuff=background_subbuff,
                              unidec_settings=unidec_settings)

        msalign_entries = self.parse_topfd_msalign()

        replacements = self.find_unidec_replacements(
            msalign_entries,
            mass_ppm_tolerance=mass_ppm_tolerance,
            charge_range=match_charge_range,
            mass_threshold=mass_threshold,
            rt_replication_window_min=rt_replication_window_min,
            min_rt_occurrences=min_rt_occurrences,
        )

        output_file = self.write_corrected_msalign(replacements, msalign_entries)

        print("\n" + "=" * 60)
        print("Pipeline complete!")
        print(f"Output: {output_file}")
        print("=" * 60)
        return output_file


def run_batch(input_dir,
              output_root,
              **pipeline_kwargs):
    """
    Process every mzML / ms2.msalign pair found in input_dir.

    For each mzML file the function looks for a matching msalign using the
    TopFD default naming convention:
        <stem>.mzML  ->  <stem>_ms2.msalign   (preferred)
        <stem>.mzML  ->  <stem>.msalign        (fallback)

    Each file pair is processed independently.  Outputs are written to a
    per-file subdirectory under output_root:
        output_root/<stem>/unidec_corrected_ms2.msalign

    Files whose output already exists are skipped, so an interrupted batch
    can be resumed by running the same command again.

    Args:
        input_dir:        Directory containing mzML and msalign files
        output_root:      Root directory for all outputs
        pipeline_kwargs:  Passed straight through to run_full_pipeline()
                          (e.g. mass_threshold=..., unidec_settings=...)
    """
    input_dir   = Path(input_dir)
    output_root = Path(output_root)

    mzml_files = sorted(set(input_dir.glob("*.mzML")) | set(input_dir.glob("*.mzml")))
    if not mzml_files:
        print(f"No mzML files found in {input_dir}")
        return

    print(f"Found {len(mzml_files)} mzML file(s) in {input_dir}")

    succeeded, failed = [], []

    for mzml_path in mzml_files:
        stem = mzml_path.stem  # e.g. "12_2D_HEPES_25ug_cut3"

        # Locate matching msalign
        msalign_path = input_dir / f"{stem}_ms2.msalign"
        if not msalign_path.exists():
            msalign_path = input_dir / f"{stem}.msalign"
        if not msalign_path.exists():
            print(f"\nWARNING: No msalign found for {mzml_path.name} — skipping")
            failed.append(mzml_path.name)
            continue

        file_output_dir = output_root / stem

        # Skip if already completed
        if (file_output_dir / OUTPUT_FILENAME).exists():
            print(f"  Skipping {mzml_path.name} — output already exists")
            succeeded.append(mzml_path.name)
            continue

        file_output_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"Processing: {mzml_path.name}")
        print(f"  msalign : {msalign_path.name}")
        print(f"  output  : {file_output_dir}")

        try:
            processor = UniDecDDAProcessor(mzml_path, file_output_dir, msalign_path)
            processor.run_full_pipeline(**pipeline_kwargs)
            succeeded.append(mzml_path.name)
        except Exception as e:
            print(f"  ERROR processing {mzml_path.name}: {e}")
            failed.append(mzml_path.name)

    print(f"\n{'='*60}")
    print(f"Batch complete: {len(succeeded)} succeeded, {len(failed)} failed")
    if failed:
        print("  Failed files:")
        for f in failed:
            print(f"    {f}")


def build_parser() -> argparse.ArgumentParser:
    d = DEFAULTS
    u = DEFAULT_UNIDEC_SETTINGS

    p = argparse.ArgumentParser(
        description="Replace TopFD precursor assignments with UniDec-derived "
                    "MS1 masses and write a corrected msalign for TopPIC.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    single = p.add_argument_group("single-file mode")
    single.add_argument("--mzml", help="Input mzML file")
    single.add_argument("--msalign", help="TopFD *_ms2.msalign file")
    single.add_argument("--output-dir", help="Output directory")

    batch = p.add_argument_group("batch mode")
    batch.add_argument("--input-dir",
                       help="Directory of <stem>.mzML and <stem>_ms2.msalign pairs")
    batch.add_argument("--output-root",
                       help="Root output directory (one subfolder per file)")

    m = p.add_argument_group("matching and filtering")
    m.add_argument("--mass-threshold", type=float, default=d['mass_threshold'],
                   help="Minimum UniDec mass (Da) to trigger replacement")
    m.add_argument("--ppm-tolerance", type=float, default=d['mass_ppm_tolerance'],
                   help="PPM tolerance for midpoint matching and replication")
    m.add_argument("--match-charge-range", type=int, nargs=2, metavar=("MIN", "MAX"),
                   default=d['match_charge_range'],
                   help="Charge states tested when matching the window midpoint")
    m.add_argument("--rt-window-min", type=float, default=d['rt_replication_window_min'],
                   help="RT window (min) for replication corroboration")
    m.add_argument("--min-rt-occurrences", type=int, default=d['min_rt_occurrences'],
                   help="Occurrences (including itself) required")

    ud = p.add_argument_group("UniDec deconvolution")
    ud.add_argument("--unidec-mass-range", type=float, nargs=2, metavar=("MIN", "MAX"),
                    default=d['unidec_mass_range'], help="Mass axis range (Da)")
    ud.add_argument("--unidec-charge-range", type=int, nargs=2, metavar=("MIN", "MAX"),
                    default=d['unidec_charge_range'], help="Charge range")
    ud.add_argument("--peak-thresh", type=float, default=d['peak_thresh'],
                    help="Peak picker threshold (fraction of tallest peak)")
    ud.add_argument("--peak-window", type=float, default=d['peak_window'],
                    help="Peak picker minimum mass separation (Da)")
    ud.add_argument("--max-peaks-per-scan", type=int, default=d['max_peaks_per_scan'],
                    help="Reject MS1 scans with more picked peaks than this")
    ud.add_argument("--background-subtype", type=int, default=d['background_subtype'],
                    help="UniDec background subtraction (0 off, 1 line, 2 curved)")
    ud.add_argument("--background-subbuff", type=float, default=d['background_subbuff'],
                    help="UniDec background subtraction parameter")
    ud.add_argument("--mz-range", type=float, nargs=2, metavar=("MIN", "MAX"),
                    default=(u['minmz'], u['maxmz']), help="m/z axis range")
    ud.add_argument("--mz-bins", type=float, default=u['mzbins'], help="m/z bin size")
    ud.add_argument("--mass-bins", type=float, default=u['massbins'],
                    help="Mass bin size (Da)")
    ud.add_argument("--psig", type=float, default=u['psig'], help="UniDec psig")
    ud.add_argument("--beta", type=float, default=u['beta'], help="UniDec beta")
    ud.add_argument("--mzsig", type=float, default=u['mzsig'], help="UniDec mzsig")

    p.add_argument("--n-workers", type=int, default=d['n_workers'],
                   help="Parallel worker processes (default: CPU count - 1)")
    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    single_mode = args.mzml or args.msalign or args.output_dir
    batch_mode  = args.input_dir or args.output_root
    if single_mode and batch_mode:
        parser.error("Use either single-file options or batch options, not both.")
    if single_mode and not (args.mzml and args.msalign and args.output_dir):
        parser.error("Single-file mode needs --mzml, --msalign and --output-dir.")
    if batch_mode and not (args.input_dir and args.output_root):
        parser.error("Batch mode needs --input-dir and --output-root.")
    if not (single_mode or batch_mode):
        parser.error("Provide single-file options or batch options.")

    pipeline_kwargs = dict(
        mass_threshold=args.mass_threshold,
        mass_ppm_tolerance=args.ppm_tolerance,
        match_charge_range=tuple(args.match_charge_range),
        unidec_mass_range=tuple(args.unidec_mass_range),
        unidec_charge_range=tuple(args.unidec_charge_range),
        peak_thresh=args.peak_thresh,
        peak_window=args.peak_window,
        max_peaks_per_scan=args.max_peaks_per_scan,
        n_workers=args.n_workers,
        rt_replication_window_min=args.rt_window_min,
        min_rt_occurrences=args.min_rt_occurrences,
        background_subtype=args.background_subtype,
        background_subbuff=args.background_subbuff,
        unidec_settings=dict(
            minmz=args.mz_range[0],
            maxmz=args.mz_range[1],
            mzbins=args.mz_bins,
            massbins=args.mass_bins,
            psig=args.psig,
            beta=args.beta,
            mzsig=args.mzsig,
        ),
    )

    if single_mode:
        processor = UniDecDDAProcessor(args.mzml, args.output_dir, args.msalign)
        processor.run_full_pipeline(**pipeline_kwargs)
    else:
        run_batch(args.input_dir, args.output_root, **pipeline_kwargs)


if __name__ == "__main__":
    main()