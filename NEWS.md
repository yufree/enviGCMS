# enviGCMS 0.5.0

- CRAN

## Major changes from last CRAN version 0.3.4

- use list to communicate results between function instead of xcms objects
- getdoe function to handle data for various experimental design
- add functions and shiny application for Short-Chain Chlorinated Paraffins analysis
- update the visulization function
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
- reformate code by formatR package

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
- add rlaplot and ridagesplot to show the distribution of the data

## Minor changes

- revise citation keywords

# enviGCMS 0.4.0

## Major changes

- Add getdoe and remove former DoE related function
- Change the plot function from xcmsset based to list based
- Add support for xcms 3 function

## Tips

- Former related functions could be found in xsetplus package here: https://github.com/yufree/xsetplus
