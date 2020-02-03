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
