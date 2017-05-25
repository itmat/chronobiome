#! /bin/bash

#Usage: path/to/chronobiome.sh "path/to/chronobiome/directory"

#Read scripts directory specified by user (if given)
SCRIPTS_DIR="$1"
if [ $# -eq 0 ]; then
        SCRIPTS_DIR="."
fi

#Date for R checkpoint package
CHECKPOINT_DATE="2016-12-20"

# Install necessary dependencies
Rscript $SCRIPTS_DIR/install.R $CHECKPOINT_DATE

# Process/combine activity,  communication, BP/heartrate, and activity datasets
Rscript $SCRIPTS_DIR/Process_Activity_Communication_BP.R $CHECKPOINT_DATE

# Generate variance explained plots and scatterplots of each variable against 
# each other variable
Rscript $SCRIPTS_DIR/Variance_Correlation_Plot.R $CHECKPOINT_DATE

# Plot principle components for activity/communication, metabolite, microbe, 
# and protein data
Rscript $SCRIPTS_DIR/PCAPlot.R $CHECKPOINT_DATE
