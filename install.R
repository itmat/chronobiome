# install.R
# Amy Campbell 2016

# Installs dependencies for the chronobiome  analysis pipeline
# Install and load checkpoint package to keep track of package versioning

library('methods')

mirror <- "http://cran.us.r-project.org"

install.packages("checkpoint", repos=mirror)

library("checkpoint")

dir.create('.checkpoint')
checkpoint("2016-12-20", checkpointLocation='.')

library("readr")
library("ggplot2")
library("Hmisc")
library("reshape2")
library("dplyr")
library("chron")
library("stats")
library("readr")
library("grid")
library("gridExtra")

sink("session.info.txt")
sessionInfo()
sink()
