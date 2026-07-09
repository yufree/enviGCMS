# enviGCMS 0.9.0

## Breaking changes

- Remove all deprecated functions from v0.7.x (32 functions including `getdata`, `getdata2`, `getmzrt2`, `svabatch`, `svaplot`, etc.). Users should migrate to current API.
- Remove `Mode()` (statistical mode function) -- unused internally since `getmdg()` moved to pmd package.
- Standardize API naming conventions:
  - Intensity column renamed from `ins` to `intensity` in `getMSP()` and `writeMSP()` spectra data.frames.
  - Rename `Integration()` to `integration()`, `GetIntegration()` to `getintegration()`, `Getisotopologues()` to `getisotopologues()` for consistent lowercase naming.

## New features (merged from MSxplorer)

- Add `HRMF()` for high-resolution mass filtering with forward/reverse scoring and Figure of Merit (FoM), based on Kwiecien et al. (2015). Uses Rdisop instead of rcdk/rJava -- no new hard dependencies.
- Add `getHRMF()` for batch HRMF processing of entire MSP files.
- Add `plotEIC()` for extracted ion chromatogram plotting from mzML/mzXML files using RaMS.
- Add `plotTopMS1Peaks()` interactive Shiny app to extract top MS1 ions from MS2 chromatograms.
- Add `plotTopMS2Peaks()` interactive Shiny app to extract top MS2 ions from MS1 chromatograms.

## Bug fixes

- Fix `getmdh()`: `cbind.data.frame(mz, MD1, MD1)` corrected to `cbind.data.frame(mz, MD1, MD2)` for ceiling method.
- Fix `runMDPlot()`: filename case mismatch (`MDPlot.rmd` -> `mdplot.Rmd`).
- Fix `getdoe()`: `sum(NROW(lv) > 1) != 0` (always TRUE) corrected to `NROW(lv) > 1`.
- Fix `getMSP()`: extract duplicated peak-parsing code into internal helper.
- Fix `getmzrtcsv()`: remove broken cross-reference to deleted `getmzrt`.
- Fix `integration()`: RT/noise window subsetting compared `> upper & < lower` (always empty), corrected to `> lower & < upper` so the function returns a real peak area.
- Fix `getmdh()`: `floor` method with more than three units assigned to `MD1_3` instead of `MD3`, causing an "object 'MD3' not found" error.
- Fix `getintegration()`: `m <= 2` branch referenced an undefined `t` and mis-parenthesised the per-scan `delta_t`; now computes `(last - first) / (n - 1)`.
- Fix `findlipid()`: `mode` branch logic was inverted -- `pos`/`neg` applied no adduct correction while the neutral branch held an unreachable `ifelse`. `pos` now corrects `[M+H]+` and `neg` corrects `[M-H]-` to neutral mass; `none` is unchanged.
- Fix `gifmr()`: default `imputation = "i"` was not a valid `getimputation()` method (silently skipped imputation); changed to `"l"`.
- Register `c.mzrt()` as an S3 method so `c()` dispatches on `mzrt` objects.
- `getcsv()`: warn instead of silently doing nothing when `type` contains none of `m`/`a`/`p`/`o`.
- Remove dead code in `getpn()` (two discarded subset expressions).

## Dependency changes

- Move `animation` from Imports to Suggests (only used by `gifmr()`).
- Move `mixtools` from Imports to Suggests (only used by `getimputation()`).
- Add `RaMS` to Suggests for lightweight mzML/mzXML access (replaces MSnbase).
- Hard Imports reduced from 11 to 8 packages.

## Code quality

- Add unit tests across 10 test files (was 1 test), including regression tests for the bug fixes above.
- Refactor `HRMF()` from 285-line monolith into 135-line function with 4 internal helpers.
- Remove 809 lines of dead deprecated code and 35 orphan man pages.

# enviGCMS 0.8.0

- remove the dependency of xcms package
- fix ppm issue in getpn

# enviGCMS 0.7.4

- CRAN
- remove rdisop related example
- update citation
- remove dependency of xcms and msnbase
- add onAttach message for bioc packages

# enviGCMS 0.7.3

- fix getalign2 to accept mz vector without ordering
- fix getalign ppm bug

# enviGCMS 0.7.2

- adope MoNA msp format for getMSP

# enviGCMS 0.7.1

- CRAN

# enviGCMS 0.7.0

- update MDplotR shiny documents
- change class to inherits
- fix save issue

# enviGCMS 0.6.9

- update writeMSP function
- add findpfc to find potential PFCs
- add getpn to connect pos and neg mode data and return mzrt object
- add support for xcms 3 EIC extraction
- add getcompare to align multiple mzrt objects

# enviGCMS 0.6.8

- goodpractice check
- spell check

# enviGCMS 0.6.7

- fix bug in getpower function
- change default par parameters for plotmrc
- remove the export of deprecated function
- fix the typo in documents
- improve getMSP to fit MoNA MSP format
- add getalign function to align two peaks vectors and remove getoverlapmass, getoverlaprt
- modify getmass to accept negative value
- add getms1anno for MS1 annotation by local compounds database with predefined paired mass distances
- add support to save target data csv without change default name
- add retention time alignment function for peaks lists
- add getalign2 function to align mass to charge ratio, retention time for the same peak list, this function could be used to refine peaks list
- fix color issues in plotrla and plotridge function

# enviGCMS 0.6.6

- CRAN

# enviGCMS 0.6.5

- change to formal group information for mzrt class
- fix log infinite issue in plotridge
- update list data

# enviGCMS 0.6.4

- plotpeak for intensity of peaks across samples or samples across peaks
- plotridge for ridgeline density plot
- plotrug for 1-d density for multiple samples

# enviGCMS 0.6.3

- add findlipid to annotate lipid class based on referenced Kendrick mass defect

# enviGCMS 0.6.2

- add getMSP to read in EI-MS or MS/MS msp file as list for annotation
- add xrankanno and dotpanno function for MS/MS X Rank and dot product annotation
- add getrangecsv for MS2 mgf file extraction

# enviGCMS 0.6.1

- remove dependency of genefilter
- remove dependency of reshape2
- remove dependency of ggplot2
- remove dependency of ggridges

# enviGCMS 0.6.0

- remove dependency of rcdk package 
- CRAN

# enviGCMS 0.5.9

- fix `...` issue

# enviGCMS 0.5.8

- add functions to get and plot density weighted intensity
- add pooled QC vignette
- add plotcc to plot calibration curve

# enviGCMS 0.5.7

- CRAN

# enviGCMS 0.5.6

- fix group issue in mzrt object
- all the group info will be imported as character while user could change the character into dataframe
- export getmzrtcsv 
- add findmet to export metabolites for certain compounds based on mass defect
- add demo data for TBBPA metabolites from this [publication](https://doi.org/10.1021/acs.est.9b02122)

# enviGCMS 0.5.5

- CRAN

# enviGCMS 0.5.4

## Major changes

- move dependances of xcms and MSnbase to suggest and remove the export for those functions
- introduce parallel computation in getdoe
- add demo data and organize examples
- add Relative Log Abundance (RLA) plots and Relative Log Abundance Ridge(RLA) plots
- add getcsv to save the list as csv file
- add getfilter to filter the list
- add getpower to deal with power analysis in metabolomics
- rewrite getmzrt and getdoe to make analysis easier
- deprecated getupload, getfeaturest and getfeaturesanova
- update vignettes
- add support for group dim larger than 2

# enviGCMS 0.5.2

## Major changes

- combine xcms and xcms 3 object function
- add support for single group plot
- add shiny application for mass defect analysis
- add function to compute high order mass defect
- add function to compare two peak list by overlap
- add function to output csv file from list object
- add function to perform pmd analysis
- add function to get chemical formula
- fix wrong object name in getmzrt/getmzrt2
- add options to export EIC object
- add function to screen organohalogen compounds
- add density plot for multiple samples

# enviGCMS 0.5.0

- CRAN

## Major changes from last CRAN version 0.3.4

- use list to communicate results between function instead of xcms objects
- getdoe function to handle data for various experimental design
- add functions and shiny application for Short-Chain Chlorinated Paraffins analysis
- update the visualization function
- support xcms 3 new objects

# enviGCMS 0.4.5

## Major changes

- remove GlobalStd function and shiny application into pmd package
- remove batch correction and simulation function into mzrtsim package
- add shiny application for sccp analysis
- add deprecated function
- update vignettes

# enviGCMS 0.4.4

## Major changes

- Shiny application for GlobalStd
- New GlobalStd algorithm
- New paired analysis based on retention time hierarchical clustering, paired mass differences(PMD), PMD frequency analysis
- Remove mass defect related analysis
- Plot function for paired analysis
- Plot function for GlobalStd result
- Organize the codes

## Minor changes

- Add Mode function
- reformat code by formatR package

# enviGCMS 0.4.3

## Major changes

- RT cluster analysis for tentative isotope, adducts, and neutral loss peaks detection
- Correlation analysis to select the feature peaks within rt groups
- Mass defect(cluster) analysis for homologous series detection
- Add csv file generation for simulation data
- Add heatmap for mzrt profile
- Add index to plotpca

## Minor changes

- Fit xcms 3 new function in getdata2
- Change the default name for metaboanalyst
- remove nonascll in CITATION
- remove faahKO package and MSnbase version for cran check

# enviGCMS 0.4.2

## Major changes

- update Vignette
- add function for SCCP analysis
- add peak list filter function based mass difference

## Minor changes

- add seealso

# enviGCMS 0.4.1

## Major changes

- add batch correction methods with p-value and q-value: sva, isva
- remove svadata and svaupload
- add mzrtsim to simulate mzrt profile
- add simmzrt to make simulation input data
- add rlaplot and ridgesplot to show the distribution of the data

## Minor changes

- revise citation keywords

# enviGCMS 0.4.0

## Major changes

- Add getdoe and remove former DoE related function
- Change the plot function from xcmsset based to list based
- Add support for xcms 3 function
