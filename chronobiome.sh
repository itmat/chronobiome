#! /bin/bash

# Install necessary dependencies
Rscript install.R

# Process/combine activity,  communication, BP/heartrate, and activity datasets
Rscript Process_Activity_Communication_BP.R

# Generate variance explained plots and scatterplots of each variable against 
# each other variable
Rscript Variance_Correlation_Plot.R

# Plot principle components for activity/communication, metabolite, microbe, 
# and protein data
Rscript PCAPlot.R
