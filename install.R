# install.R
# Amy Campbell 2016

# Installs dependencies for the chronobiome  analysis pipeline
# Install and load checkpoint package to keep track of package versioning


mirror <- "http://cran.us.r-project.org"
install.packages("methods", repos=mirror)
install.packages("checkpoint", repos=mirror)

library('methods')
library("checkpoint")


dir.create('.checkpoint')
checkpoint("2016-12-20", checkpointLocation='.')

packages <- c(
  "readr",
  "ggplot2",
  "Hmisc",
  "reshape2",
  "dplyr",
  "chron",
  "stats",
  "readr",
  "grid",
  "gridExtra"
  )
sink("session.info.txt")
sessionInfo()
sink()
