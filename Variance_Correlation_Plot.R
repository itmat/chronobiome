################################################
# Amy Campbell 2016-2017
#
# Plots variability explained by variables in activity, communication,
# energy, and blood pressure/heart rate data 
# 
# Makes scatterplots of each combination of activity, communication, 
# blood pressure/heart rate, and energy
################################################

###############
# Load Packages
###############

# Read date of snapshot to use
args = commandArgs(TRUE)
checkpoint_date = args[1]

# Use checkpoint created by install script
library("checkpoint")
checkpoint(snapshotDate=checkpoint_date, checkpointLocation='.')

library("ggplot2")
library("Hmisc")
library("reshape2")
library("stats")
library("readr")
library("chron")
library("grid")
library("gridExtra")


##################
# Define functions
##################


GetR2 <- function(var1, var2) {
  # Generates a simple linear regression model of var1~var2, and returns
  # R-squared for that model to represent proportion of V1 explained by V2
  # :param: var1 - Variable to 'be explained' by var2
  # :param: var2 - Variable to 'explain' var1
  Var1name <- colnames(var1)[1]
  Var2name <-  colnames(var2)[1]
  model <- lm(as.numeric(var1) ~ as.numeric(var2))
  return(c(Var1name, Var2name, (summary(model)$r.squared)))
}
  

PlotR2s <- function(chrono_df, title, Sequence) {
  # Generates geom_tile plot of the % variance explained (R2) for each variable
  # by each other variable.
  # Also outputs matrix of R-squared values in .csv format with the same 
  # :param: chrono_df - Dataframe to plot, including all variables(activity, 
  #                     communication, and in certain time frames,
  #                     blood pressure) to be compared in this geom_tile plot 
  #                     note: each variable in the dataframe will be compared 
  #                     to each other variable in the dataframe (and itself);
  #                     because these will be compared using R-squared values 
  #                     yielded from simple linear regression models, they 
  #                     must be numeric 
  # :param: title - title to display on plot 
  # :param: Sequence - string indicating which order of ints to use to plot 
  #                    variables on the variance explained output plot:
  #                    sequence.actcom (just activity/communication variables)
  #                    sequence.actcomBP (activity & communication
  #                                       with blood pressure/HR data)
  #                    sequence.energy (for datasets including activity,
  #                                     communication, blood pressure,
  #                                     and energy data)
  #                                     
  chrono_df["Subject"] <- NULL
  chrono_df["TimeSubjectIndex"] <- NULL            
  r2Matrix <- c()
  for (a in colnames(chrono_df)) {
    for (b in colnames(chrono_df)) {
      r2Matrix <-
        rbind(r2Matrix,
              GetR2(as.matrix(chrono_df[a]),
                    as.matrix(chrono_df[b])
              ))
    }
  }
  plot.matrix <- data.frame(r2Matrix)
  colnames(plot.matrix) <- c("factor.1", "factor.2", "variability.explained")
  plot.matrix$variability.explained <-
    as.numeric(as.character(plot.matrix$variability.explained))
  # Order variables using their alphabetically-determined level indices such 
  # that activraphy, actigraphy circadian, 
  # communication, communication circadian, and biometric variables are
  # grouped together in the output chart
  # sequence.actcomBP <- c(8, 9, 10, 7, 21, 26, 25, 24, 5, 19, 18, 17, 1,
  #                       2, 3, 6, 20, 22, 13, 16, 15, 14, 23, 4, 11, 27, 12)
  # sequence.actcom <- c(7, 8, 9 , 6, 18, 22, 21, 20, 4, 16, 15, 14, 1,
  #                      2, 3, 5, 17, 19, 10, 13, 12, 11)
  # sequence.energy <- c(8, 9, 10, 7, 23, 28, 27, 26, 5, 21, 20, 18, 1, 2, 3,
  #                     6, 22, 19, 17, 24, 13, 16, 15, 14, 25, 4, 11, 29, 12)
  # sequence.4monthenergy <- c(7, 8, 9, 6, 20, 24, 23, 22, 4, 18, 17, 15, 1, 2,
  #                           3, 5, 19, 16, 14, 21, 10, 13, 12, 11)
  
  # if (Sequence == "ActCom"){
  #     f1 <- factor(plot.matrix$factor.1,
  #                  levels(plot.matrix$factor.1)[sequence.actcom])
  #     f2 <- factor(plot.matrix$factor.2,
  #                  levels(plot.matrix$factor.2)[sequence.actcom])           
  # } else if (Sequence == "BP") {
  #     f1 <- factor(plot.matrix$factor.1,
  #                  levels(plot.matrix$factor.1)[sequence.actcomBP])
  #     f2 <- factor(plot.matrix$factor.2,
  #                  levels(plot.matrix$factor.2)[sequence.actcomBP])
  # } else if (Sequence == "energy") {
  #     f1 <- factor(plot.matrix$factor.1,
  #                  levels(plot.matrix$factor.1)[sequence.energy])
  #     f2 <- factor(plot.matrix$factor.2,
  #                  levels(plot.matrix$factor.2)[sequence.energy])
  # } else {
  #     f1 <- factor(plot.matrix$factor.1,
  #                  levels(plot.matrix$factor.1)[sequence.4monthenergy])
  #     f2 <- factor(plot.matrix$factor.2,
  #                  levels(plot.matrix$factor.2)[sequence.4monthenergy]) 
  # }
  
  # Or alternatively, order the variables by explicitly using the names of each
  # factor level, rather than the number corresponding to its alphabetical order.
  # The factor levels are still hard-coded, but this makes it a little easier
  # to incorporate new values.
  sequence.actcomBP <- c("Communication.amplitude", "Communication.period",
                         "Communication.phase", "circadian.signal.com",
                         "log.signal.com", "sqrt.Interaction.Diversity",
                         "SMS.Length", "SMS.Count", "Call.Count",
                         "MobilityRadius.amplitude", "MobilityRadius.period",
                         "MobilityRadius.phase", "circadian.signal.mobR",
                         "log.Mobility.Radius", "Mobility.amplitude", "Mobility.period",
                         "Mobility.phase", "circadian.signal.mob", "log.Mobility",
                         "Lux.amplitude", "Lux.period", "Lux.phase", "circadian.signal.lux",
                         "log.Luminosity", "activity.amplitude", "activity.period",
                         "activity.phase", "circadian.signal.act", "log.signal.act",
                         "log.Steps", "log.activity.Vector.Magnitude", "log.Axis3",
                         "log.Axis2", "log.Axis1", "pulse.pressure", "arterial.pressure",
                         "diastolic.bp", "systolic.bp", "heart.rate")
  sequence.actcom <- c("Communication.amplitude", "Communication.period",
                       "Communication.phase", "circadian.signal.com",
                       "log.signal.com", "sqrt.Interaction.Diversity",
                       "SMS.Length", "SMS.Count", "Call.Count",
                       "MobilityRadius.amplitude", "MobilityRadius.period",
                       "MobilityRadius.phase", "circadian.signal.mobR",
                       "log.Mobility.Radius", "Mobility.amplitude", "Mobility.period",
                       "Mobility.phase", "circadian.signal.mob", "log.Mobility",
                       "Lux.amplitude", "Lux.period", "Lux.phase", "circadian.signal.lux",
                       "log.Luminosity", "activity.amplitude", "activity.period",
                       "activity.phase", "circadian.signal.act", "log.signal.act",
                       "log.Steps", "log.activity.Vector.Magnitude", "log.Axis3",
                       "log.Axis2", "log.Axis1")
  sequence.energy  <- c("Communication.amplitude", "Communication.period",
                        "Communication.phase", "circadian.signal.com",
                        "log.signal.com", "sqrt.Interaction.Diversity",
                        "SMS.Length", "SMS.Count", "Call.Count",
                        "MobilityRadius.amplitude", "MobilityRadius.period",
                        "MobilityRadius.phase", "circadian.signal.mobR",
                        "log.Mobility.Radius", "Mobility.amplitude", "Mobility.period",
                        "Mobility.phase", "circadian.signal.mob", "log.Mobility",
                        "Lux.amplitude", "Lux.period", "Lux.phase", "circadian.signal.lux",
                        "log.Luminosity", "activity.amplitude", "activity.period",
                        "activity.phase", "circadian.signal.act", "log.signal.act",
                        "log.MET.rate", "log.kcals", "log.Steps",
                        "log.activity.Vector.Magnitude", "log.Axis3", "log.Axis2",
                        "log.Axis1", "pulse.pressure", "arterial.pressure",
                        "diastolic.bp", "systolic.bp", "heart.rate")
  sequence.4monthenergy <- c("Communication.amplitude", "Communication.period",
                             "Communication.phase", "circadian.signal.com",
                             "log.signal.com", "sqrt.Interaction.Diversity",
                             "SMS.Length", "SMS.Count", "Call.Count",
                             "MobilityRadius.amplitude", "MobilityRadius.period",
                             "MobilityRadius.phase", "circadian.signal.mobR",
                             "log.Mobility.Radius", "Mobility.amplitude", "Mobility.period",
                             "Mobility.phase", "circadian.signal.mob", "log.Mobility",
                             "Lux.amplitude", "Lux.period", "Lux.phase", "circadian.signal.lux",
                             "log.Luminosity", "activity.amplitude", "activity.period",
                             "activity.phase", "circadian.signal.act", "log.signal.act",
                             "log.MET.rate", "log.kcals", "log.Steps",
                             "log.activity.Vector.Magnitude", "log.Axis3", "log.Axis2",
                             "log.Axis1")
  

  if (Sequence == "ActCom"){
    f1 <- factor(plot.matrix$factor.1, sequence.actcom)
    f2 <- factor(plot.matrix$factor.2, sequence.actcom)
  } else if (Sequence == "BP") {
    f1 <- factor(plot.matrix$factor.1, sequence.actcomBP)
    f2 <- factor(plot.matrix$factor.2, sequence.actcomBP)
  } else if (Sequence == "energy") {
    f1 <- factor(plot.matrix$factor.1, sequence.energy)
    f2 <- factor(plot.matrix$factor.2, sequence.energy)
  } else {
    f1 <- factor(plot.matrix$factor.1, sequence.4monthenergy)
    f2 <- factor(plot.matrix$factor.2, sequence.4monthenergy)
  }
  # generate geom_tile plot
  plot <-
    ggplot(data = plot.matrix,
           aes(x = f1, y = f2, fill = variability.explained)) +
           geom_tile(color = "black") +
           scale_fill_gradient(low = "white", high = "steelblue4") +
           theme(axis.text.x = element_text(angle = 90, size = 8,
                                            hjust = 1, vjust = 0.5),
           axis.text.y = element_text(size = 8)) + ggtitle(title) +
           xlab("Factor 1") + ylab("Factor 2")

  return(plot)
}

ParseSubject <- function(row, half) {
  # Function to be applied to all rows of 
  # the dataframe in question
  # :param: row - the entire row of the dataframe
  # :param: half - indicates whether the "subject" information
  # is in the first or second half of the string
  timesubjectindex <- row["TimeSubjectIndex"]
  strlist <- (strsplit(toString(timesubjectindex), "_"))
  
  if (toString(half) == "2") {
    return(unlist(strlist)[2])
  } else {
    return(as.numeric(unlist(strlist)[1]))
  }
}


TimeList_window1 <- c(seq(922, 1232))
TimeList_window2 <- c(seq(1328, 1545))

TimeList <- c(TimeList_window1, TimeList_window2)
###############################
# Heatmap of variance explained
###############################

transformed.communication.activity <- 
  readr::read_csv("transformed.communication.activity.csv")
transformed.communication.activity["X1"] <- NULL

# Two 48 Hour visits containing measurements of blood pressure, heart rate
BP_HR <- readr::read_csv("heartrate.bp.csv")
BP_HR["X1"] <- NULL
act.com.bp.hr.vars4months <- dplyr::full_join(BP_HR, transformed.communication.activity, 
                                              by = "TimeSubjectIndex")
bp.HR.com.act <- na.omit(act.com.bp.hr.vars4months)
postscript("Variability.Act.Com.BP.eps")
PlotR2s(bp.HR.com.act[,5:dim(bp.HR.com.act)[2]], 
        paste0("Variance Explained in Activity,", 
              "Communication, Biometric Data (Visits 1 and 2)"), "BP")
dev.off()

# 48 hour visit 1
Visit1_bp.HR.com.act <- subset(bp.HR.com.act, Times %in% TimeList_window1)

postscript("Varability.Act.Com.BP.Visit1.Feb.eps")
PlotR2s(Visit1_bp.HR.com.act[, 5:dim(Visit1_bp.HR.com.act)[2]],
        " Variance Explained in Activity, Communication, Biometric Data(Visit 1)",
        "BP")
dev.off()

# 48 hour visit 2
Visit2_bp.HR.com.act <- subset(bp.HR.com.act, Times %in% TimeList_window2)
postscript("Varability.Act.Com.BP.Visit2.Feb.eps")
PlotR2s(Visit2_bp.HR.com.act[, 5:dim(Visit2_bp.HR.com.act)[2]],
        "Variance Explained in Activity, Communication, Biometric Data (Visit 2)", "BP")
dev.off()

energy <- read.csv("energy.csv")
energy["X"] <- NULL
energy["Times"] <- apply(energy, 1, ParseSubject, 1)

Energy.4months <- dplyr::full_join(transformed.communication.activity, energy,
                                   by = "TimeSubjectIndex")
Energy.4months <- na.omit(Energy.4months)

Full.with.energy <- dplyr::full_join(bp.HR.com.act, energy,
                                     by = "TimeSubjectIndex")
Full.with.energy <- na.omit(Full.with.energy)

postscript("Variance_Explained_WithEnergy.eps")
# Plot R^2s excluding Times variable
PlotR2s(Full.with.energy[, 6:dim(Full.with.energy)[2] - 2],
        paste0("Variance Explained in Activity, Communication, Blood Pressure,",
              " and Energy Variables"),
        # "4monthenergy") #I think this should be "energy" instead of "4monthenergy"
        "energy")         #since Full.with.energy includes the BP stats.
dev.off()

# All four months by subject (no blood pressure/heart rate/ ActCom data)
transformed.communication.activity["Subject"] <- 
  apply(transformed.communication.activity, 1, ParseSubject, 2)

# Subject HCR001 
HCR001.4months <- subset(transformed.communication.activity, Subject == "HCR001")
postscript("HCR001_Variance.4months.eps",  width = 480, height = 480)
PlotR2s(HCR001.4months, "Variance Explained in HCR001 (all 4 months)",
        "ActCom")
dev.off()

# HCR003
HCR003.4months <- subset(transformed.communication.activity, Subject == "HCR003")
postscript("HCR003_Variance.4months.eps",  width = 480, height = 480)
PlotR2s(HCR003.4months, "Variance Explained in HCR003 (all 4 months)",
        "ActCom")
dev.off()

# HCR004
HCR004.4months <- subset(transformed.communication.activity, Subject == "HCR004") 
postscript("HCR004_Variance.4months.eps",  width = 480, height = 480)
PlotR2s(HCR004.4months, "Variance Explained in HCR004 (all 4 months)",
        "ActCom")
dev.off()

# HCR006
HCR006.4months <- subset(transformed.communication.activity, Subject == "HCR006")
postscript("HCR006_Variance.4months.eps",  width = 480, height = 480)
PlotR2s(HCR006.4months, "Variance Explained in HCR006 (all 4 months)",
        "ActCom")
dev.off()

# HCR008
HCR008.4months <- subset(transformed.communication.activity, Subject == "HCR008")
postscript("HCR008_Variance.4months.eps",  width = 480, height = 480)
PlotR2s(HCR008.4months, "Variance Explained in HCR008 (all 4 months)",
        "ActCom")
dev.off()

# HCR009
HCR009.4months <- subset(transformed.communication.activity, Subject == "HCR009")
postscript("HCR009_Variance.4months.eps",  width = 480, height = 480)
PlotR2s(HCR009.4months, "Variance Explained in HCR009 (all 4 months)",
        "ActCom")
dev.off()


# All four months (no blood pressure/heart rate data)
postscript("Variability.Activity.Communication.eps",  width = 480, height = 480)
PlotR2s(transformed.communication.activity,
        "Variance Explained in Activity and Communication Data", "ActCom")
dev.off()


######################################

# All four months by subject (with energy)
Energy.4months["Subject"] <- (apply(Energy.4months, 1, ParseSubject, 2))
Energy.4months["X1"] <- NULL
Energy.4months["Times"] <- NULL
# Subject HCR001 
HCR001.4months <- subset(Energy.4months, Subject == "HCR001")
postscript("HCR001_Variance.4months_withEnergy.eps",
           width = 480, height = 480)
PlotR2s(HCR001.4months, "Variance Explained in HCR001 (all 4 months)",
        "4monthenergy")
dev.off()

# HCR003
HCR003.4months <- subset(Energy.4months, Subject == "HCR003")
postscript("HCR003_Variance.4months_withEnergy.eps",
           width = 480, height = 480)
PlotR2s(HCR003.4months, "Variance Explained in HCR003 (all 4 months)",
        "4monthenergy")
dev.off()

#HCR003
HCR004.4months <- subset(Energy.4months, Subject == "HCR004")
postscript("HCR004_Variance.4months_withEnergy_withEnergy.eps",
           width = 480, height = 480)
PlotR2s(HCR004.4months, "Variance Explained in HCR004 (all 4 months)",
        "4monthenergy")
dev.off()

# HCR006
HCR006.4months <- subset(Energy.4months, Subject == "HCR006")
postscript("HCR006_Variance.4months_withEnergy.eps", width = 480, height = 480)
PlotR2s(HCR006.4months, "Variance Explained in HCR006 (all 4 months)",
        "4monthenergy")
dev.off()

# HCR008
HCR008.4months <- subset(Energy.4months, Subject == "HCR008")
postscript("HCR008_Variance.4months_withEnergy.eps", width = 480, height = 480)
PlotR2s(HCR008.4months, "Variance Explained in HCR008 (all 4 months)",
        "4monthenergy")
dev.off()

# HCR009
HCR009.4months <- subset(Energy.4months, Subject == "HCR009")
postscript("HCR009_Variance.4months_withEnergy.eps", width = 480, height = 480)
PlotR2s(HCR009.4months, "Variance Explained in HCR009 (all 4 months)",
        "4monthenergy")
dev.off()

postscript("Variability.Activity.Communication_WithEnergy.eps",
           width = 480, height = 480)
PlotR2s(Energy.4months,
        "Variance Explained in Activity and Communication Data (4 months)",
        "4monthenergy")
dev.off()

# Both Visits by Subject all variables 
Full.with.energy["Subject"] <- (apply(Full.with.energy, 1, ParseSubject, 2))

# Subject HCR001 
HCR001.set <- subset(Full.with.energy, Subject == "HCR001")
postscript("HCR001_Variance_BothVisits.eps",  width = 480, height = 480)
PlotR2s(HCR001.set[, 6:dim(HCR001.set)[2] - 2], "Variance Explained in HCR001",
        "energy")
dev.off()

# HCR003
HCR003.set <- subset(Full.with.energy, Subject == "HCR003")
postscript("HCR003_Variance_BothVisits.eps",  width = 480, height = 480)
PlotR2s(HCR003.set[, 6:dim(HCR003.set)[2] - 2], "Variance Explained in HCR003",
        "energy")
dev.off()

# HCR004
HCR004.set <- subset(Full.with.energy, Subject == "HCR004")
postscript("HCR004_Variance_BothVisits.eps",  width = 480, height = 480)
PlotR2s(HCR004.set[, 6:dim(HCR004.set)[2] - 2], "Variance Explained in HCR004",
        "energy")
dev.off()

# HCR006
HCR006.set <- subset(Full.with.energy, Subject == "HCR006")
HCR006.set["Subject"] <- NULL
postscript("HCR006_Variance_BothVisits.eps",  width = 480, height = 480)
PlotR2s(HCR006.set[, 6:dim(HCR006.set)[2] - 2], "Variance Explained in HCR006",
        "energy")
dev.off()
  
# HCR008
HCR008.set <- subset(Full.with.energy, Subject == "HCR008")
postscript("HCR008_Variance_BothVisits.eps",  width = 480, height = 480)
PlotR2s(HCR008.set[, 6:dim(HCR008.set)[2] - 2], "Variance Explained in HCR008",
        "energy")
dev.off()

# HCR009
HCR009.set <- subset(Full.with.energy, Subject == "HCR009")
postscript("HCR009_Variance_BothVisits.eps",  width = 480, height = 480)
PlotR2s(HCR009.set[, 6:dim(HCR009.set)[2] - 2], "Variance Explained in HCR009",
        "energy")
dev.off()


################################
# VISIT 1
# Subject HCR001 
HCR001.set1 <- subset(subset(Full.with.energy, Subject == "HCR001"), Days <= 43)
postscript("HCR001_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR001.set1[, 6:dim(HCR001.set1)[2] - 2],
        "Variance Explained in HCR001 (visit 1)", "energy")
dev.off()

# HCR003
HCR003.set1 <- subset(subset(Full.with.energy, Subject == "HCR003"), (Days <= 43))
postscript("HCR003_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR003.set1[, 6:dim(HCR003.set1)[2] - 2],
        "Variance Explained in HCR003 (visit 1)", "energy")
dev.off()

# HCR004
HCR004.set1 <- subset(subset(Full.with.energy, Subject == "HCR004"), (Days<=45))
postscript("HCR004_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR004.set1[, 6:dim(HCR004.set1)[2] - 2],
        "Variance Explained in HCR004 (visit 1)", "energy")
dev.off()

# HCR006
HCR006.set1 <- subset(subset(Full.with.energy, Subject == "HCR006"), (Days<=44))
postscript("HCR006_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR006.set1[, 6:dim(HCR006.set1)[2] - 2],
        "Variance Explained in HCR006 (visit 1)", "energy")
dev.off()

# HCR008
HCR008.set1 <- subset(subset(Full.with.energy, Subject == "HCR008"), Days<=45)
postscript("HCR008_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR008.set1[, 6:dim(HCR008.set1)[2] - 2],
        "Variance Explained in HCR008 (visit 1)", "energy")
dev.off()

# HCR009
HCR009.set1 <- subset(subset(Full.with.energy, Subject == "HCR009"), Days <= 51)
postscript("HCR009_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR009.set1[, 6:dim(HCR009.set1)[2] - 2],
        "Variance Explained in HCR009 (visit 1)", "energy")
dev.off()


# VISIT 2
# Subject HCR001 
HCR001.set2 <- subset(subset(Full.with.energy, Subject == "HCR001"), Days>=55)
postscript("HCR001_Variance_visit2.eps",  width = 480, height = 480)
PlotR2s(HCR001.set2[, 6:dim(HCR001.set2)[2] - 2],
        "Variance Explained in HCR001 (visit 2)", "energy")
dev.off()

# HCR003
HCR003.set2 <- subset(subset(Full.with.energy, Subject == "HCR003"), Days >= 55)
postscript("HCR003_Variance_visit2.eps",  width = 480, height = 480)
PlotR2s(HCR003.set2[, 6:dim(HCR003.set2)[2] - 2],
        "Variance Explained in HCR003 (visit 2)", "energy")
dev.off()

# HCR004
HCR004.set2 <- subset(subset(Full.with.energy, Subject == "HCR004"), Days >= 55)
postscript("HCR004_Variance_visit2.eps",  width = 480, height = 480)
PlotR2s(HCR004.set2[, 6:dim(HCR004.set2)[2] - 2],
        "Variance Explained in HCR004 (visit 2)", "energy")
dev.off()

# HCR006
HCR006.set2 <- subset(subset(Full.with.energy, Subject == "HCR006"), Days >= 57)
postscript("HCR006_Variance_visit2.eps",  width = 480, height = 480)
PlotR2s(HCR006.set2[, 6:dim(HCR006.set2)[2] - 2],
        "Variance Explained in HCR006 (visit 2)", "energy")
dev.off()

# HCR008
HCR008.set2 <- subset(subset(Full.with.energy, Subject == "HCR008"), Days >= 55)
# HCR008.set2 is empty for the full set of variables including energy as a result of
# missing values from visit 2.

# HCR009
HCR009.set2 <- subset(subset(Full.with.energy, Subject == "HCR009"), Days >= 62)
postscript("HCR009_Variance_visit2.eps",  width = 480, height = 480)
PlotR2s(HCR009.set2[, 6:dim(HCR009.set2)[2] - 2],
        "Variance Explained in HCR009 (visit 2)", "energy")
dev.off()

##############
# Scatterplots 
##############
# Use Full.with.energy to save scatterplots of every combination of variables
# within the 48 hour window (406 plots)

# Save subject info
subj <- Full.with.energy["Subject"]
Full.with.energy.reduced <- subset(Full.with.energy,
                                   select = -c(Days, Times.y,
                                               TimeSubjectIndex, Times.x, Subject))

# Get a list of all 29-choose-2 combinations of two variables
combos <-  gtools::combinations(n = 29, r = 2,
                                v = colnames(Full.with.energy.reduced),
                                repeats.allowed = F)
Full.with.energy.reduced["Subject"] <- subj
# Save indices in the 'combos' list from 1 to 406
sequence <- seq(1, dim(combos)[1])

# Break up combos list into groups of 4 (and one group of 2)
# Save to "scatterplots.pdf" with 4 of the 406 plots are on
# each page 
chunks <- split(sequence, ceiling(seq_along(sequence) / 4))

colorblind_Palette <- c("#000000", "#0072B2", "#56B4E9", "#F0E442",
                        "#D55E00", "#CC79A7")

pdf("scatterplots.pdf")
for (chunk in chunks) {
  if (length(chunk) == 4) {
    pair1 <- combos[chunk[1], ]
    pair2 <- combos[chunk[2], ]
    pair3 <- combos[chunk[3], ]
    pair4 <- combos[chunk[4], ]
    
    plot1 <- qplot(x = Full.with.energy.reduced[pair1[1]],
                   y = Full.with.energy.reduced[pair1[2]],
                   data = Full.with.energy.reduced, xlab = pair1[1], ylab = pair1[2],
                   colour = Subject) + scale_colour_manual(values=colorblind_Palette)
    plot2 <- qplot(x = Full.with.energy.reduced[pair2[1]],
                   y = Full.with.energy.reduced[pair2[2]],
                   data = Full.with.energy.reduced, xlab = pair2[1], ylab = pair2[2], 
                   colour = Subject) + scale_colour_manual(values=colorblind_Palette)
    plot3 <- qplot(x = Full.with.energy.reduced[pair3[1]],
                   y = Full.with.energy.reduced[pair3[2]],
                   data = Full.with.energy.reduced, xlab = pair3[1], ylab = pair3[2],
                   colour = Subject) + scale_colour_manual(values=colorblind_Palette)
    plot4 <- qplot(x = Full.with.energy.reduced[pair4[1]],
                   y = Full.with.energy.reduced[pair4[2]],
                   data = Full.with.energy.reduced, xlab = pair4[1], ylab = pair4[2], 
                   colour = Subject) + scale_colour_manual(values=colorblind_Palette)
    grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2)

  } else {
    plot1 <- qplot(x = Full.with.energy.reduced[pair1[1]],
                   y = Full.with.energy.reduced[pair1[2]],
                   data = Full.with.energy.reduced, xlab = pair1[1], ylab = pair1[2],
                   colour = Subject) + scale_colour_manual(values=colorblind_Palette)
    plot2 <- qplot(x = Full.with.energy.reduced[pair2[1]],
                   y = Full.with.energy.reduced[pair2[2]],
                   data = Full.with.energy.reduced, xlab = pair2[1], ylab = pair2[2],
                   colour = Subject) + scale_colour_manual(values=colorblind_Palette)
    grid.arrange(plot1, plot2, ncol = 2, nrow = 2)
  }
}

dev.off()

