Chronobiome Data Analyses
==============================================
**Amy Campbell, 2016**
---------------------------------------------------------------------------------
Repo Contents
---------------------------------------------------------------------------------
* **chronobiome.sh**
  - Shell script, runs the scripts below
* **install.R**
   - Installs necessary R dependencies
   - Saves record of installation with [checkpoint](https://cran.r-project.org/web/packages/checkpoint/index.html) package
   - Outputs:
  - *session.info.txt*
* **Process_Activity_Communication_BP.R** 
   - Assigns common hour/subject index to all datapoints in activity, communication, blood pressure/HR, and energy data
   - Transforms activity and communication variables to improve normality
   - Outputs:
 - *transformed.community.activity.csv*
 - *community.activity.csv* 
 - *heartrate.bp.csv*
 - *energy.csv*
* **PCAPlot.R**
   - Performs PCA on combined activity and communication data across all 4 months
   - Plots activity and communication PCA, and plots metabolomic and microbiomic PCAs performed by Seth Rhodes
   - Outputs: 
 - *Act_Com_PCAPlotAll.eps*
 - *Act_Com_PCAPlot_48.eps* 
 - *Facets_Act_Com_PC2PC3.eps*
 - *Facets_Act_Com_PC1PC3.eps*
 - *Facets_Act_Com_PC2PC3.eps*
 - *Saliva_Metabs_PCA.eps*
 - *Saliva_Microbes_PCA.eps*
 - *Plasma_Proteins_PCA.eps*
 - *Plasma_Metabs_PCA.eps*
 - *Macronutrients.eps*
 - *Variance_Explained_WithEnergy.eps*
* **Variance_Correlation_Plot.R**
  - Generates R-squared (variability explained) heatmaps for activity & communication variables for all 4 months and  for activity, communication, and blood pressure / heart rate variables
  - Generates .pdf of scatterplots for each activity, communication, energy, and blood pressure/heartrate variable
against each other 
  - Outputs:
    - *HCR001_Variance.4months.eps*
    - *HCR003_Variance.4months.eps*
    - *HCR004_Variance.4months.eps*
    - *HCR006_Variance.4months.eps*
    - *HCR008_Variance.4months.eps*
    - *HCR009_Variance.4months.eps*
    - *Variability.Activity.Communication.eps*
    - *HCR001_Variance.4months_withEnergy*
    - *HCR003_Variance.4months_withEnergy*
    - *HCR004_Variance.4months_withEnergy*
    - *HCR006_Variance.4months_withEnergy*
    - *HCR008_Variance.4months_withEnergy*
    - *HCR009_Variance.4months_withEnergy*
    - *Variability.Activity.Communication_WithEnergy.eps*
    - *HCR001_Variance_BothVisits*
    - *HCR003_Variance_BothVisits*
    - *HCR004_Variance_BothVisits*
    - *HCR006_Variance_BothVisits*
    - *HCR008_Variance_BothVisits*
    - *HCR009_Variance_BothVisits*
    - *HCR001_Variance_visit1*
    - *HCR003_Variance_visit1*
    - *HCR004_Variance_visit1*
    - *HCR006_Variance_visit1*
    - *HCR008_Variance_visit1*
    - *HCR009_Variance_visit1*
    - *HCR001_Variance_visit2*
    - *HCR003_Variance_visit2*
    - *HCR004_Variance_visit2*
    - *HCR006_Variance_visit2*
    - *HCR009_Variance_visit2*
 
