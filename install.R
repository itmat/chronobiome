# install.R
# Amy Campbell 2016

# Installs dependencies for the chronobiome  analysis pipeline

# Install and load checkpoint package to keep track of package versioning
install.packages("methods")
install.packages("checkpoint", repos = "http://cran.us.r-project.org")

library("checkpoint")
dir.create(".checkpoint")
checkpoint("2016-12-20")

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
install.packages(packages)
