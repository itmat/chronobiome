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
* **Plot.R**
   - Performs PCA on combined activity and communication data across all 4 months
   - Generates R-squared (variability explained) heatmaps for activity & communication variables for all 4 months and for activity, communication, and blood pressure / heart rate variables 
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
 
