#### UniDec Precursor Correction for TopFD / TopPIC

Replaces TopFD precursor assignments with masses derived from MS1 deconvolution by Unidec, and writes a corrected `msalign` file that can be searched with TopPIC.

\[citation]

## 

For large proteins, isotopic resolution is often poor or absent, so features may be missed or mis-assigned by TopFD, and the resulting precursor mass is wrong or missing. This script deconvolutes the parent MS1 scan with UniDec (which does not require isotopic resolution), finds the UniDec mass that matches the isolated mass, and rewrites the precursor fields in the TopFD `msalign` entry. The fragment ions are left untouched.

#### 

##### Workflow

```
raw  ->  mzML (centroid)
              |
              +--> TopFD --> <stem>\_ms2.msalign ---+
              |                                    |
              +------------------------------------+--> this script --> unidec\_corrected\_ms2.msalign --> TopPIC
```

##### Requirements

Tested with the following versions.

|Component|Version|
|-|-|
|Python|3.11.5|
|pyOpenMS|3.5.0|
|NumPy|1.26.4|
|UniDec (Python package)|8.1|
|TopFD|1.8|
|TopPIC|1.8|
|mzML conversion|MSConvert, Proteowizard 3.0|
|OS|Windows 11|

## 

##### Usage

###### Single file

```bash
python unidec\_precursor\_correction.py \\
    --mzml sample.mzML \\
    --msalign sample\_ms2.msalign \\
    --output-dir results/sample
```

###### Batch (all mzML/msalign pairs in a directory)

```bash
python unidec\_precursor\_correction.py \\
    --input-dir data/ \\
    --output-root results/ \\
    --n-workers 8
```

Batch mode writes one subfolder per file and skips files whose output already exists, so an interrupted run can be resumed by repeating the command.

#### Parameters

The defaults are the values used in the publication. Everything below can be set on the command line (`--help` lists all flags) or as keyword arguments to `run\_full\_pipeline()` / `run\_batch()`.

###### Matching and filtering

|CLI flag|Python argument|Default|Description|
|-|-|-|-|
|`--mass-threshold`|`mass\_threshold`|`30000`|Minimum UniDec mass (Da) for a replacement. Smaller masses are left to TopFD.|
|`--ppm-tolerance`|`mass\_ppm\_tolerance`|`100`|PPM tolerance when matching the window midpoint to a UniDec mass. Also used for the replication check.|
|`--match-charge-range`|`match\_charge\_range`|`5 115`|Charge states tried when converting the window midpoint m/z to a mass.|
|`--rt-window-min`|`rt\_replication\_window\_min`|`1.0`|RT window (minutes) for the replication check.|
|`--min-rt-occurrences`|`min\_rt\_occurrences`|`2`|Number of MS2 spectra (including itself) that must show the same mass within the RT window.|

###### UniDec deconvolution

|CLI flag|Python argument|Default|Description|
|-|-|-|-|
|`--unidec-mass-range`|`unidec\_mass\_range`|`4000 110000`|Mass axis range (Da).|
|`--unidec-charge-range`|`unidec\_charge\_range`|`5 115`|Charge range searched by UniDec.|
|`--peak-thresh`|`peak\_thresh`|`0.1`|Peak picker threshold, as a fraction of the tallest peak in the mass spectrum.|
|`--peak-window`|`peak\_window`|`15.0`|Minimum mass separation (Da) between distinct picked peaks.|
|`--max-peaks-per-scan`|`max\_peaks\_per\_scan`|`15`|Scans yielding more picked peaks than this are rejected as too noisy or unresolved, and none of their masses are used.|
|`--background-subtype`|`background\_subtype`|`2`|UniDec background subtraction: 0 off, 1 line, 2 curved.|
|`--background-subbuff`|`background\_subbuff`|`100`|Background subtraction parameter.|
|`--mz-range`|`unidec\_settings\['minmz'/'maxmz']`|`500 3000`|m/z range used by UniDec.|
|`--mz-bins`|`unidec\_settings\['mzbins']`|`1.0`|m/z bin size.|
|`--mass-bins`|`unidec\_settings\['massbins']`|`1.0`|Mass bin size (Da).|
|`--psig`|`unidec\_settings\['psig']`|`10.0`|UniDec `psig`.|
|`--beta`|`unidec\_settings\['beta']`|`100`|UniDec `beta`.|
|`--mzsig`|`unidec\_settings\['mzsig']`|`0.5`|UniDec `mzsig`.|
|*(not exposed)*|n/a|`psfun = 0`|Peak shape function, fixed to Gaussian. Editable in `DEFAULT\_UNIDEC\_SETTINGS`.|

###### Other

|CLI flag|Python argument|Default|Description|
|-|-|-|-|
|`--n-workers`|`n\_workers`|CPU count - 1|Parallel worker processes for deconvolution.|



#### Output

For each input, the script writes `unidec\_corrected\_ms2.msalign` into the output directory.

The file contains only the spectra that were replaced, not the whole TopFD output. For each of these, the five precursor fields are rewritten and everything else is copied from the original entry:

|Field|New value|
|-|-|
|`PRECURSOR\_MZ`|Isolation window midpoint|
|`PRECURSOR\_CHARGE`|Charge at which the window midpoint matched the UniDec mass|
|`PRECURSOR\_MASS`|UniDec mass|
|`PRECURSOR\_INTENSITY`|UniDec peak height|
|`PRECURSOR\_FEATURE\_ID`|`-1` (no TopFD feature)|

## 

#### Citing

If you use this software, please cite:

* \[citation]
* UniDec: Marty, M. T. *et al.* Bayesian deconvolution of mass and ion mobility spectra: from binary interactions to polydisperse ensembles. *Anal. Chem.* 87, 4370–4376 (2015).
* TopFD / TopPIC: \[add the TopFD and TopPIC references you used]
* pyOpenMS: Röst, H. L. *et al.* OpenMS: a flexible open-source software platform for mass spectrometry data analysis. *Nat. Methods* 13, 741–748 (2016).

