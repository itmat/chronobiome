################################################
# Amy Campbell 2016-2017
#
# Fits and plots principal component analysis of 
# activity and communication variables
#
# Plots variability explained by variables in 
# activity, communication, and blood pressure/
# heart rate data 
################################################

###############
# Load Packages
###############

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


GetR2 <- function(var1, var2){
  # Generates a simple linear regression model of var1~var2, and returns
  # R-squared for that model to represent proportion of V1 explained by V2
  # :param: var1 - Variable to 'be explained' by var2
  # :param: var2 - Variable to 'explain' var1
  Var1name <- colnames(var1)[1]
  Var2name <-  colnames(var2)[1]
  model <- lm(as.numeric(var1) ~ as.numeric(var2))
  return(c(Var1name, Var2name, (summary(model)$r.squared)))
}
  

PlotR2s <- function(chrono_df, title, Sequence){
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
  #                    sequence_actcom (just activity/communication variables)
  #                    sequence_actcomBP (activity & communication
  #                                       with blood pressure/HR data)
  #                    sequence_energy (for datasets including activity,
  #                                     communication, blood pressure,
  #                                     and energy data)
  #                                     
              
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
  
  PlotMatrix <- data.frame(r2Matrix)
  colnames(PlotMatrix) <- c("factor.1", "factor.2", "variability.explained")
  PlotMatrix$variability.explained <- as.numeric(as.character(
    PlotMatrix$variability.explained))
  # Order variables using their alphabetically-determined level indices such 
  # that activraphy, actigraphy circadian, 
  # communication, communication circadian, and biometric variables are
  # grouped together in the output chart
  sequence_actcomBP = c(8, 9, 10, 7, 21, 26, 25, 24, 5, 19, 18, 17, 1,
                        2, 3, 6, 20, 22, 13, 16, 15, 14, 23, 4, 11, 27, 12)
  sequence_actcom = c(7, 8, 9 , 6, 18, 22, 21, 20, 4, 16, 15, 14, 1,
                       2, 3, 5, 17, 19, 10, 13, 12, 11)
  sequence_energy = c(8, 9, 10, 7, 23, 28, 27, 26, 5, 21, 20, 18, 1, 2, 3,
                      6, 22, 19, 17, 24, 13, 16, 15, 14, 25, 4, 11, 29, 12)
  sequence_4monthenergy = c(7, 8, 9, 6, 20, 24, 23, 22, 4, 18, 17, 15, 1, 2,
                            3, 5, 19, 16, 14, 21, 10, 13, 12, 11)

  if (Sequence == "ActCom"){
    f1 = factor(PlotMatrix$factor.1,
                levels(PlotMatrix$factor.1)[sequence_actcom])
    f2 = factor(PlotMatrix$factor.2,
                levels(PlotMatrix$factor.2)[sequence_actcom])           
  }
  # Order variables for 48 hour communication/activity/BP dataset
  else if (Sequence == "BP"){
    f1 <- factor(PlotMatrix$factor.1,
                 levels(PlotMatrix$factor.1)[sequence_actcomBP])
    f2 <- factor(PlotMatrix$factor.2,
                 levels(PlotMatrix$factor.2)[sequence_actcomBP])
  }
  else if (Sequence == "energy"){
    f1 <- factor(PlotMatrix$factor.1,
                 levels(PlotMatrix$factor.1)[sequence_energy])
    f2 <- factor(PlotMatrix$factor.2,
                 levels(PlotMatrix$factor.2)[sequence_energy])
  } else {
    f1 <- factor(PlotMatrix$factor.1,
                 levels(PlotMatrix$factor.1)[sequence_4monthenergy])
    f2 <- factor(PlotMatrix$factor.2,
                 levels(PlotMatrix$factor.2)[sequence_4monthenergy]) 
  }
  # generate geom_tile plot
  plot <-
    ggplot(data = PlotMatrix,
           aes(x = f1, y = f2, fill = variability.explained)) +
           geom_tile(color = "black") +
           scale_fill_gradient(low="white", high = "steelblue4") +
           theme(axis.text.x = element_text(angle = 90, size = 8) +
           axis.text.y = element_text(size = 8)) + ggtitle(title) +
           xlab("Factor 1") + ylab("Factor 2")

  return(plot)
}

parsesubject <- function(row, half){
  # Function to be applied to all rows of 
  # the dataframe in question
  # :param: row - the entire row of the dataframe
  # :param: half - indicates whether the "subject" information
  # is in the first or second half of the string
  timesubjectindex <- row["TimeSubjectIndex"]
  strlist <- (strsplit(toString(timesubjectindex), "_"))
  
  if(toString(half) == "2"){
    return(unlist(strlist)[2])
  }
  else{
    return(as.numeric(unlist(strlist)[1]))
  }
}

TimeList_window1 <- c(seq(922, 1232))
TimeList_window2 <- c(seq(1328, 1545))

###############################
# Heatmap of variance explained
###############################

transformed.communication.activity <- 
  readr::read_csv("transformed.communication.activity.csv")
transformed.communication.activity["X1"] = NULL

# Two 48 Hour visits containing measurements of blood pressure, heart rate
BP_HR <- readr::read_csv("heartrate.bp.csv")
Full <- dplyr::full_join(BP_HR, transformed.communication.activity,
                         by = "TimeSubjectIndex")
bp.HR.com.act <- na.omit(Full)
postscript("Variability.Act.Com.BP.eps")
PlotR2s(bp.HR.com.act[,5:dim(bp.HR.com.act)[2]], 
        paste("Variance Explained in Activity,", 
              " Communication, Biometric Data (Visits 1 and 2)", "BP"))

dev.off()


# 48 hour visit 1
Visit1_bp.HR.com.act <- subset(bp.HR.com.act, Times %in% TimeList_window1)

postscript("Varability.Act.Com.BP.Visit1.Feb.eps")
PlotR2s(Visit1_bp.HR.com.act[, 5:dim(Visit1_bp.HR.com.act)[2]],
        "Variance Explained in Activity, Communication, Biometric Data(Visit 1)",
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
energy["Times"] <- apply(energy, 1, parsesubject, 1)

Energy.4months <- dplyr::full_join(transformed.communication.activity, energy,
                                   by = "TimeSubjectIndex")
Energy.4months <- na.omit(Energy.4months)

Full.with.energy <- dplyr::full_join(bp.HR.com.act, energy,
                                     by = "TimeSubjectIndex")
Full.with.energy <- na.omit(Full.with.energy)
Full.with.energy["Subject"] <- NULL

postscript("Variance_Explained_WithEnergy.eps")
PlotR2s(Full.with.energy[, 6:dim(Full.with.energy)[2] - 1],
        paste("Variance Explained in Activity, Communication, Blood Pressure,",
              " and Energy Variables"),
        "somethingelse")
dev.off()

# All four months (no blood pressure/heart rate data)
transformed.communication.activity["TimeSubjectIndex"] <- NULL
transformed.communication.activity["Subject"] <- NULL
postscript("Variability.Activity.Communication.eps",  width = 480, height = 480)
PlotR2s(transformed.communication.activity,
        "Variance Explained in Activity and Communication Data", "ActCom")
dev.off()


######################################

# All four months by subject (with energy)
Energy.4months["Subject"] <- (apply(Energy.4months, 1, parsesubject, 2))
Energy.4months["X1"] <- NULL
Energy.4months["Times"] <- NULL
# Subject HCR001 
HCR001_4months <- subset(Energy.4months, Subject == "HCR001")
HCR001_4months["Subject"] <- NULL
HCR001_4months["TimeSubjectIndex"] <- NULL
postscript("HCR001_Variance_4months_withEnergy.eps",
           width = 480, height = 480)
PlotR2s(HCR001_4months, "Variance Explained in HCR001 (all 4 months)",
        "4monthenergy")
dev.off()

# HCR003
HCR003_4months <- subset(Energy.4months, Subject=="HCR003")
HCR003_4months["Subject"] <- NULL
HCR003_4months["TimeSubjectIndex"] <- NULL
postscript("HCR003_Variance_4months_withEnergy.eps",
           width = 480, height = 480)
PlotR2s(HCR003_4months, "Variance Explained in HCR003 (all 4 months)",
        "4monthenergy")
dev.off()

#HCR003
HCR004_4months <- subset(Energy.4months, Subject == "HCR004")
HCR004_4months["Subject"] <- NULL
HCR004_4months["TimeSubjectIndex"] <- NULL

postscript("HCR004_Variance_4months_withEnergy_withEnergy.eps",
           width = 480, height = 480)
PlotR2s(HCR004_4months, "Variance Explained in HCR004 (all 4 months)",
        "4monthenergy")
dev.off()

# HCR006
HCR006_4months <- subset(Energy.4months, Subject == "HCR006")
HCR006_4months["Subject"] <- NULL
HCR006_4months["TimeSubjectIndex"] <- NULL

postscript("HCR006_Variance_4months_withEnergy.eps", width = 480, height = 480)
PlotR2s(HCR006_4months, "Variance Explained in HCR006 (all 4 months)",
        "4monthenergy")
dev.off()

# HCR008
HCR008_4months <- subset(Energy.4months, Subject=="HCR008")
HCR008_4months["Subject"] <- NULL
HCR008_4months["TimeSubjectIndex"] <- NULL
postscript("HCR008_Variance_4months_withEnergy.eps", width = 480, height = 480)
PlotR2s(HCR008_4months, "Variance Explained in HCR008 (all 4 months)",
        "4monthenergy")
dev.off()

# HCR009
HCR009_4months <- subset(Energy.4months, Subject=="HCR009")
HCR009_4months["Subject"] <- NULL
HCR009_4months["TimeSubjectIndex"] <- NULL

postscript("HCR009_Variance_4months_withEnergy.eps", width = 480, height = 480)
PlotR2s(HCR009_4months, "Variance Explained in HCR009 (all 4 months)",
        "4monthenergy")
dev.off()

Energy.4months["Subject"] <- NULL
Energy.4months["TimeSubjectIndex"] <- NULL
postscript("Variability.Activity.Communication_WithEnergy.eps",
           width = 480, height = 480)
PlotR2s(Energy.4months,
        "Variance Explained in Activity and Communication Data (4 months)",
        "4monthenergy")
dev.off()
View(Energy.4months)
####


# All four months by subject (no blood pressure/heart rate/ ActCom data)
transformed.communication.activity["Subject"] <- 
  (apply(transformed.communication.activity, 1, parsesubject, 2))
# Subject HCR001 
HCR001_4months <- subset(transformed.communication.activity, Subject == "HCR001")
HCR001_4months["Subject"] <- NULL
HCR001_4months["TimeSubjectIndex"] <- NULL
postscript("HCR001_Variance_4months.eps",  width = 480, height = 480)
PlotR2s(HCR001_4months, "Variance Explained in HCR001 (all 4 months)",
        "ActCom")
dev.off()

# HCR003
HCR003_4months <- subset(transformed.communication.activity, Subject == "HCR003")
HCR003_4months["Subject"] <- NULL
HCR003_4months["TimeSubjectIndex"] <- NULL
postscript("HCR003_Variance_4months.eps",  width = 480, height = 480)
PlotR2s(HCR003_4months, "Variance Explained in HCR003 (all 4 months)",
        "ActCom")
dev.off()

#HCR003
HCR004_4months <- subset(transformed.communication.activity, Subject == "HCR004")
HCR004_4months["Subject"] <- NULL
HCR004_4months["TimeSubjectIndex"] <- NULL

postscript("HCR004_Variance_4months.eps",  width = 480, height = 480)
PlotR2s(HCR004_4months, "Variance Explained in HCR004 (all 4 months)",
        "ActCom")
dev.off()

# HCR006
HCR006_4months <- subset(transformed.communication.activity, Subject == "HCR006")
HCR006_4months["Subject"] <- NULL
HCR006_4months["TimeSubjectIndex"] <- NULL

postscript("HCR006_Variance_4months.eps",  width = 480, height = 480)
PlotR2s(HCR006_4months, "Variance Explained in HCR006 (all 4 months)",
        "ActCom")
dev.off()

# HCR008
HCR008_4months <- subset(transformed.communication.activity, Subject == "HCR008")
HCR008_4months["Subject"] <- NULL
HCR008_4months["TimeSubjectIndex"] <- NULL
postscript("HCR008_Variance_4months.eps",  width = 480, height = 480)
PlotR2s(HCR008_4months, "Variance Explained in HCR008 (all 4 months)",
        "ActCom")
dev.off()

# HCR009
HCR009_4months <- subset(transformed.communication.activity, Subject=="HCR009")
HCR009_4months["Subject"] <- NULL
HCR009_4months["TimeSubjectIndex"] <- NULL

postscript("HCR009_Variance_4months.eps",  width = 480, height = 480)
PlotR2s(HCR009_4months, "Variance Explained in HCR009 (all 4 months)",
        "ActCom")
dev.off()


# Both Visits by Subject all variables 
Full.with.energy["Subject"] <- (apply(Full.with.energy, 1, parsesubject, 2))

# Subject HCR001 
HCR001_set <- subset(Full.with.energy, Subject=="HCR001")
HCR001_set["Subject"] <- NULL
postscript("HCR001_Variance_BothVisits.eps",  width = 480, height = 480)
PlotR2s(HCR001_set[, 6:dim(HCR001_set)[2]-1], "Variance Explained in HCR001",
        "energy")
dev.off()

# HCR003
HCR003_set <- subset(Full.with.energy, Subject == "HCR003")
HCR003_set["Subject"] <- NULL
postscript("HCR003_Variance_BothVisits.eps",  width = 480, height = 480)
PlotR2s(HCR003_set[, 6:dim(HCR003_set)[2] - 1], "Variance Explained in HCR003",
        "energy")
dev.off()
  
#HCR004
HCR004_set <- subset(Full.with.energy, Subject == "HCR004")
HCR004_set["Subject"] <- NULL
postscript("HCR004_Variance_BothVisits.eps",  width = 480, height = 480)
PlotR2s(HCR004_set[, 6:dim(HCR004_set)[2] - 1], "Variance Explained in HCR004",
        "energy")
dev.off()

# HCR006
HCR006_set <- subset(Full.with.energy, Subject == "HCR006")
HCR006_set["Subject"] <- NULL
postscript("HCR006_Variance_BothVisits.eps",  width = 480, height = 480)
PlotR2s(HCR006_set[, 6:dim(HCR006_set)[2] - 1], "Variance Explained in HCR006",
        "energy")
dev.off()
  
# HCR008
HCR008_set <- subset(Full.with.energy, Subject == "HCR008")
HCR008_set["Subject"] <- NULL
postscript("HCR008_Variance_BothVisits.eps",  width = 480, height = 480)
PlotR2s(HCR008_set[, 6:dim(HCR008_set)[2]-1], "Variance Explained in HCR008",
        "energy")
dev.off()


# HCR009
HCR009_set <- subset(Full.with.energy, Subject == "HCR009")
HCR009_set["Subject"] <- NULL
postscript("HCR009_Variance_BothVisits.eps",  width = 480, height = 480)
PlotR2s(HCR009_set[, 5:(dim(HCR009_set)[2]-1)], "Variance Explained in HCR009",
        "energy")
dev.off()


################################
# VISIT 1
# Subject HCR001 
HCR001_set1 <- subset(subset(Full.with.energy, Subject == "HCR001"), Days<=43)
HCR001_set1["Subject"] <- NULL
postscript("HCR001_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR001_set1[, 5:(dim(HCR001_set1)[2] - 1)],
        "Variance Explained in HCR001 (visit 1)", "energy")
dev.off()

# HCR003
HCR003_set1 <- subset(subset(Full.with.energy, Subject == "HCR003"), (Days <= 43))
HCR003_set1["Subject"] <- NULL
postscript("HCR003_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR003_set1[, 5:(dim(HCR003_set1)[2] - 1)],
        "Variance Explained in HCR003 (visit 1)", "energy")
dev.off()

#HCR004
HCR004_set1 <- subset(subset(Full.with.energy, Subject == "HCR004"), (Days<=45))
HCR004_set1["Subject"] <- NULL
postscript("HCR004_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR004_set1[, 5:(dim(HCR004_set1)[2] - 1)],
        "Variance Explained in HCR004 (visit 1)", "energy")
dev.off()

# HCR006
HCR006_set1 <- subset(subset(Full.with.energy, Subject == "HCR006"), (Days<=44))
HCR006_set1["Subject"] <- NULL
postscript("HCR006_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR006_set1[, 5:(dim(HCR006_set1)[2] - 1)],
        "Variance Explained in HCR006 (visit 1)", "energy")
dev.off()

# HCR008
HCR008_set1 <- subset(subset(Full.with.energy, Subject == "HCR008"), Days<=45)
HCR008_set1["Subject"] <- NULL
postscript("HCR008_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR008_set1[, 5:(dim(HCR008_set1)[2]-1)],
        "Variance Explained in HCR008 (visit 1)", "energy")
dev.off()

# HCR009
HCR009_set1 <- subset(subset(Full.with.energy, Subject == "HCR009"), Days <= 51)
HCR009_set1["Subject"] <- NULL
postscript("HCR009_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR009_set1[, 5:(dim(HCR009_set1)[2]-1)],
        "Variance Explained in HCR009 (visit 1)", "energy")
dev.off()


# VISIT 2
# Subject HCR001 
HCR001_set2 <- subset(subset(Full.with.energy, Subject=="HCR001"), Days>=55)
HCR001_set2["Subject"] <- NULL
postscript("HCR001_Variance_visit2.eps",  width = 480, height = 480)
PlotR2s(HCR001_set2[, 5:(dim(HCR001_set2)[2]-1)],
        "Variance Explained in HCR001 (visit 2)", "energy")
dev.off()

# HCR003
HCR003_set2 <- subset(subset(Full.with.energy, Subject == "HCR003"), Days >= 55)
HCR003_set2["Subject"] <- NULL
postscript("HCR003_Variance_visit2.eps",  width = 480, height = 480)
PlotR2s(HCR003_set2[, 5:(dim(HCR003_set2)[2] - 1)],
        "Variance Explained in HCR003 (visit 2)", "energy")
dev.off()

#HCR004
HCR004_set2 <- subset(subset(Full.with.energy, Subject == "HCR004"), Days >= 55)
HCR004_set2["Subject"] <- NULL
postscript("HCR004_Variance_visit2.eps",  width = 480, height = 480)
PlotR2s(HCR004_set2[, 5:(dim(HCR004_set2)[2] - 1)],
        "Variance Explained in HCR004 (visit 2)", "energy")
dev.off()

# HCR006
HCR006_set2 <- subset(subset(Full.with.energy, Subject == "HCR006"), Days >= 57)
HCR006_set2["Subject"] <- NULL
postscript("HCR006_Variance_visit2.eps",  width = 480, height = 480)
PlotR2s(HCR006_set2[, 5:(dim(HCR006_set2)[2]-1)],
        "Variance Explained in HCR006 (visit 2)", "energy")
dev.off()

# HCR008
HCR008_set2 <- subset(subset(Full.with.energy, Subject == "HCR008"), Days >=55)
HCR008_set2["Subject"] <- NULL
postscript("HCR008_Variance_visit2.eps",  width = 480, height = 480)
PlotR2s(HCR008_set2[, 5:(dim(HCR008_set2)[2] - 1)],
        "Variance Explained in HCR008 (visit 2)", "energy")
dev.off()

# HCR009
HCR009_set2 <- subset(subset(Full.with.energy, Subject == "HCR009"), Days >= 62)
HCR009_set2["Subject"] <- NULL
postscript("HCR009_Variance_visit2.eps",  width = 480, height = 480)
PlotR2s(HCR009_set2[, 5:(dim(HCR009_set2)[2]-1)],
        "Variance Explained in HCR009 (visit 2)", "energy")
dev.off()




################################
# Communication and activity PCA
################################

communication.activity <- (read.csv("communication.activity.csv"))
communication.activity['subject'] <- 
  (apply(communication.activity, 1, parsesubject, 2))

TimeSubj <- data.frame(communication.activity$TimeSubjectIndex)
subj <- data.frame(communication.activity$subject)
length <- as.numeric(dim(communication.activity)[2])
Matrix_For_PCA <- communication.activity[,3:(length-1)]
Act_Com_PCA <- prcomp(Matrix_For_PCA, scale.= TRUE)
ActCom_PCA_Matrix <- data.frame(Act_Com_PCA$x)

ActCom_PCA_Matrix["Subject"] <- subj
write.table(Act_Com_PCA$rotation, "Act_Com_Loadings.csv", sep = ",")


ActCom_PCA_Matrix["TimeSubjectIndex"] <- TimeSubj
  
ActCom_PCA_Matrix["TimesOnly"] <- 
  apply(ActCom_PCA_Matrix, 1, parsesubject, 1)

windowed_ActCom_PCA <- 
  subset(ActCom_PCA_Matrix, TimesOnly %in% TimeList)

###########
# PCA Plots
###########
colorblind_Palette <- c("#000000", "#0072B2", "#56B4E9", "#F0E442",
                        "#D55E00", "#CC79A7")

postscript("Act_Com_PCAPlotAll.eps", width=480, height=480)
p1 <- qplot(x = PC2, y = PC3, data = ActCom_PCA_Matrix, colour = Subject) + 
            scale_colour_manual(values = colorblind_Palette)
p2 <- qplot(x = PC1, y = PC3, data = ActCom_PCA_Matrix, colour = Subject) + 
            scale_colour_manual(values = colorblind_Palette)
p3 <- qplot(x = PC1, y = PC2, data = ActCom_PCA_Matrix, colour = Subject) + 
            scale_colour_manual(values = colorblind_Palette)
grid.arrange(p1, p2, p3, ncol = 3, 
             top = paste("Principal Components for Combined Activity and ",
                         "Communication (All Measurements)"))
dev.off()

postscript("Act_Com_PCAPlot_48.eps")
p_1 <- qplot(x = PC2, y = PC3, data = windowed_ActCom_PCA, colour = Subject) +
             scale_colour_manual(values=colorblind_Palette) 
p_2 <- qplot(x = PC1, y = PC3, data = windowed_ActCom_PCA, colour = Subject) +
             scale_colour_manual(values=colorblind_Palette)
p_3 <- qplot(x = PC1, y = PC2, data = windowed_ActCom_PCA, colour = Subject) +
             scale_colour_manual(values=colorblind_Palette)
grid.arrange(p_1, p_2, p_3, ncol = 3,
             top = paste("Principal Components for Combined Activity",
                         "and Communication (Visits 1 and 2)"))
dev.off()

postscript("Facets_Act_Com_PC2PC3.eps")
fp_1 <- p_1 + facet_grid(. ~Subject)
grid.arrange(fp_1, ncol = 1,
             top = paste("Principal Components for Combined Activity",
                         "and Communication (Visits 1 and 2)"))
dev.off()

postscript("Facets_Act_Com_PC1PC3.eps")
fp_2 <- p_2 + facet_grid(. ~Subject) 
grid.arrange(fp_2, ncol = 1, 
             top = paste("Principal Components for Combined Activity and",
                         " Communication (Visits 1 and 2)"))
dev.off()

postscript("Facets_Act_Com_PC1PC2.eps")
fp_3 <- p_3 + facet_grid(. ~Subject) 
grid.arrange(fp_3, ncol = 1, 
             top = paste("Principal Components for Combined Activity and", 
                         " Communication (Visits 1 and 2)"))
dev.off()


# Load and plot other PCAs (Generated by S. Rhodes)
Saliva_Metabolites <- read.csv(
  paste("chronobiome/Seth_PCA/PC_Scores/SalivaMetabolites_BothVisits",
        "_ScoresFor10PCs.csv"))
colnames(Saliva_Metabolites) <- c("Visit", "Subject", "Hour", "PC1", "PC2",
                                  "PC3", "PC4", "PC5", "PC6", "PC7",
                                  "PC8", "PC9", "PC10")
postscript("Saliva_Metabs_PCA.eps")
sm1 <- qplot(x = PC2, y = PC3, data = Saliva_Metabolites, colour = Subject)
             + scale_colour_manual(values = colorblind_Palette)
sm2 <- qplot(x = PC1, y = PC3, data = Saliva_Metabolites, colour = Subject)
             + scale_colour_manual(values=colorblind_Palette)
sm3 <- qplot(x = PC1, y = PC2, data = Saliva_Metabolites, colour = Subject)
             + scale_colour_manual(values=colorblind_Palette)
grid.arrange(sm1, sm2, sm3, ncol = 3,
             top =
               "Principal Components for Saliva Metabolites(Visits 1 and 2)")
dev.off()

Saliva_Microbes <-
  read.csv(paste("chronobiome/Seth_PCA/PC_Scores/SalivaMicrobes_",
                 "BothVisits_ScoresFor5PCs.csv"))
colnames(Saliva_Microbes) <-(c("X", "Subject", "Hour","PC1",
                               "PC2", "PC3", "PC4", "PC5"))
postscript("Saliva_Microbes_PCA.eps")
smic1 <- qplot(x = PC2, y = PC3, data = Saliva_Microbes, colour = Subject) +
              scale_colour_manual(values = colorblind_Palette) 
smic2 <- qplot(x = PC1, y = PC3, data = Saliva_Microbes, colour=Subject) +
              scale_colour_manual(values = colorblind_Palette) 
smic3 <- qplot(x = PC1, y = PC2, data = Saliva_Microbes, colour = Subject) +
              scale_colour_manual(values = colorblind_Palette)
grid.arrange(smic1, smic2, smic3, ncol = 3,
             top = "Principal Components for Saliva Microbes(Visits 1 and 2)")
dev.off()

Plasma_Proteins <- read.csv(
  paste("chronobiome/Seth_PCA/PC_Scores/PlasmaProteins_1stVisits_HCR009-",
        "12hrRemove_ScoresFor10PCs.csv"))
Plasma_Proteins <- na.omit(Plasma_Proteins)
colnames(Plasma_Proteins) <- c("X", "Subject", "Hour", "PC1", "PC2", "PC3",
                               "PC4", "PC5", "PC6", "PC7","PC8","PC9","PC10")
postscript("Plasma_Proteins_PCA.eps")
pp1 <- qplot(x = PC2, y = PC3, data = Plasma_Proteins, colour = Subject) +
            scale_colour_manual(values = colorblind_Palette) 
pp2 <- qplot(x = PC1, y = PC3, data = Plasma_Proteins, colour = Subject) +
            scale_colour_manual(values = colorblind_Palette) 
pp3 <- qplot(x = PC1, y = PC2, data = Plasma_Proteins, colour = Subject) +
            scale_colour_manual(values = colorblind_Palette)
grid.arrange(pp1, pp2, pp3, ncol = 3,
             top = "Principal Components for Plasma Proteins(Visit 1)")
dev.off()

Plasma_Metabolites <- read.csv(
  paste("chronobiome/Seth_PCA/PC_Scores/", 
        "PlasmaMetabolites_BothVisits_ScoresFor10PCs.csv"))
colnames(Plasma_Metabolites) <- c("Visit", "Subject", "Hour", "PC1",
                                  "PC2", "PC3", "PC4", "PC5", "PC6",
                                  "PC7", "PC8", "PC9", "PC10")
postscript("Plasma_Metabs_PCA.eps")
pm1 <- qplot(x = PC2, y = PC3, data = Plasma_Metabolites, colour = Subject) +
            scale_colour_manual(values = colorblind_Palette)
pm2 <- qplot(x = PC1, y = PC3, data = Plasma_Metabolites, colour = Subject) +
            scale_colour_manual(values = colorblind_Palette) 
pm3 <- qplot(x = PC1, y = PC2, data = Plasma_Metabolites, colour = Subject) +
          scale_colour_manual(values = colorblind_Palette)
grid.arrange(pm1, pm2, pm3, ncol = 3,
             top = "Principal Components for Plasma Metabolites(Visits 1 and 2)")
dev.off()

Macronutrients <- read.csv(
  "/home/amycampbell/Downloads/MacroNutrientData_SmallerChunks_BothVisits.csv")
postscript("Macronutrients.eps")
n1 <- qplot(x = Protein, y = TotalFat, data = Macronutrients, colour = Subject) +
            scale_colour_manual(values = colorblind_Palette)
n2 <- qplot(x = Protein, y = Carbohydrate, data = Macronutrients, colour = Subject) +
            scale_colour_manual(values=colorblind_Palette)
n3 <- qplot(x = Carbohydrate, y = TotalFat, data = Macronutrients, colour = Subject) +
            scale_colour_manual(values=colorblind_Palette)
grid.arrange(n1, n2, n3, ncol = 3, top = "Macronutrients(Visits 1 and 2)")
dev.off()

Metabolites_Microbiome_BothVisits<- 
  read.csv(paste("chronobiome/Seth_PCA/PC_Scores/ScoresJan17/",
                 "MetabsAndMicrobes_BothVisits_ScoresFor5PCs.csv"))
colnames(Metabolites_Microbiome_BothVisits) <-
  c("X", "Subject", "Hour", "PC1", "PC2", "PC3", "PC4", "PC5")
postscript("Metabolites_Microbes_BothVisits.eps")
mm1 <- qplot(x = PC2, y = PC3, data = Metabolites_Microbiome_BothVisits,
             colour = Subject) +
  scale_colour_manual(values = colorblind_Palette)
mm2 <- qplot(x = PC1, y = PC3, data = Metabolites_Microbiome_BothVisits,
             colour = Subject) +
  scale_colour_manual(values=colorblind_Palette) 
mm3 <- qplot(x = PC1, y = PC2, data = Metabolites_Microbiome_BothVisits,
             colour = Subject) +
  scale_colour_manual(values = colorblind_Palette)
grid.arrange(mm1, mm2, mm3, ncol = 3,
             top = paste("Principal Components for Metabolites and Microbes",
                  "Combined (Visits 1 and 2)"))
dev.off()

Metabolites_Microbiome_Visit1 <- 
  read.csv(paste("chronobiome/Seth_PCA/PC_Scores/ScoresJan17/",
                 "MetabsAndMicrobes_Visit1_ScoresFor5PCs.csv"))
colnames(Metabolites_Microbiome_Visit1) <-
  c("X", "Subject", "Hour", "PC1", "PC2", "PC3", "PC4", "PC5")
postscript("Metabolites_Microbes_Visit1.eps")
mm11 <- qplot(x = PC2, y = PC3, data = Metabolites_Microbiome_Visit1,
              colour=Subject) + scale_colour_manual(values = colorblind_Palette)
mm12 <- qplot(x = PC1, y = PC3, data = Metabolites_Microbiome_Visit1,
              colour = Subject) + scale_colour_manual(values = colorblind_Palette) 
mm13 <- qplot(x = PC1, y = PC2, data = Metabolites_Microbiome_Visit1,
              colour = Subject) + scale_colour_manual(values = colorblind_Palette)
grid.arrange(mm11, mm12, mm13, ncol = 3,
             top = paste("Principal Components for Metabolites and Microbes",
                         "Combined (Visit 1)"))
dev.off()

Metabolites_Microbiome_Visit2 <- read.csv(
  paste("chronobiome/Seth_PCA/PC_Scores/ScoresJan17/MetabsAndMicrobes_Visit2",
        "_ScoresFor5PCs.csv"))
colnames(Metabolites_Microbiome_Visit2) <- c("X", "Subject", "Hour", "PC1",
                                             "PC2", "PC3", "PC4", "PC5")
postscript("Metabolites_Microbes_Visit2.eps")
mm21 <- qplot(x = PC2, y = PC3, data = Metabolites_Microbiome_Visit2,
              colour = Subject) +
  scale_colour_manual(values=colorblind_Palette)
mm22 <- qplot(x=PC1, y=PC3, data=Metabolites_Microbiome_Visit2,
              colour=Subject) +
  scale_colour_manual(values = colorblind_Palette) 
mm23 <- qplot(x = PC1, y=PC2, data = Metabolites_Microbiome_Visit2,
              colour = Subject) +
  scale_colour_manual(values = colorblind_Palette)
grid.arrange(mm21, mm22, mm23, ncol = 3,
             top = paste("Principal Components for Metabolites",
             " and Microbes Combined (Visit 2)"))
dev.off()

TriOmic_Visit1 <- read.csv(paste("chronobiome/Seth_PCA/PC_Scores/ScoresJan17",
                                 "/TriOmic_Visit1_ScoresFor5PCs.csv"))
colnames(TriOmic_Visit1) <- c("X", "Subject", "Hour", "PC1", "PC2",
                              "PC3", "PC4", "PC5")
postscript("TriOmic_Visit1.eps")
t1 <- qplot(x = PC2, y = PC3, data = TriOmic_Visit1, colour = Subject) +
            scale_colour_manual(values = colorblind_Palette)
t2 <- qplot(x = PC1, y = PC3, data = TriOmic_Visit1, colour = Subject) +
            scale_colour_manual(values = colorblind_Palette) 
t3 <- qplot(x = PC1, y = PC2, data = TriOmic_Visit1, colour = Subject) +
            scale_colour_manual(values = colorblind_Palette)
grid.arrange(t1, t2, t3, ncol = 3,
             top = paste("Principal Components for Metabolites,",
                         " Proteins and Microbes Combined (Visit 1)"))
dev.off()

# Plot Seth's updated PCA
Saliva_Metabolites_Jan <- 
  read.csv(paste("chronobiome/Seth_PCA/PC_Scores/ScoresJan17/",
                 "SalivaMetabolites_BothVisits_ScoresFor5PCs.csv"))
colnames(Saliva_Metabolites_Jan) <- c("Visit", "Subject", "Hour",
                                      "PC1", "PC2", "PC3", "PC4", "PC5")
postscript("Saliva_Metabs_PCA_Jan17.eps")
sm1 <- qplot(x = PC2, y = PC3, data = Saliva_Metabolites_Jan, colour = Subject) +
             scale_colour_manual(values = colorblind_Palette)
sm2 <- qplot(x = PC1, y = PC3, data = Saliva_Metabolites_Jan, colour = Subject) +
             scale_colour_manual(values = colorblind_Palette)
sm3 <- qplot(x = PC1, y = PC2, data = Saliva_Metabolites_Jan, colour = Subject) +
             scale_colour_manual(values = colorblind_Palette)
grid.arrange(sm1, sm2, sm3, ncol = 3,
             top = "Principal Components for Saliva Metabolites(Visits 1 and 2)")
dev.off()

Saliva_Microbes_Jan <- 
  read.csv(paste("chronobiome/Seth_PCA/PC_Scores/ScoresJan17/SalivaMicrobes_",
                 "BothVisits_ScoresFor5PCs.csv"))
colnames(Saliva_Microbes_Jan) <- 
  (c("X", "Subject", "Hour", "PC1", "PC2", "PC3", "PC4", "PC5"))
postscript("Saliva_Microbes_PCA_Jan.eps")
smic1 <- qplot(x = PC2, y = PC3, data = Saliva_Microbes_Jan, colour = Subject) +
              scale_colour_manual(values = colorblind_Palette) 
smic2 <- qplot(x = PC1, y = PC3, data = Saliva_Microbes_Jan, colour = Subject) +
              scale_colour_manual(values = colorblind_Palette) 
smic3 <- qplot(x = PC1, y = PC2, data = Saliva_Microbes_Jan, colour = Subject) +
              scale_colour_manual(values = colorblind_Palette)
grid.arrange(smic1, smic2, smic3, ncol = 3,
             top = "Principal Components for Saliva Microbes(Visits 1 and 2)")
dev.off()

Plasma_Proteins_Jan <- 
  read.csv(paste("chronobiome/Seth_PCA/PC_Scores/ScoresJan17/",
           "PlasmaProteins_Visit1_ScoresFor5PCs.csv"))
Plasma_Proteins_Jan <- na.omit(Plasma_Proteins_Jan)
colnames(Plasma_Proteins_Jan) <- c("X", "Subject", "Hour",
                                   "PC1", "PC2", "PC3", "PC4", "PC5")
postscript("Plasma_Proteins_PCA_Jan.eps")
pp1 <- qplot(x = PC2, y = PC3, data = Plasma_Proteins_Jan, colour = Subject) +
            scale_colour_manual(values = colorblind_Palette) 
pp2 <- qplot(x = PC1, y = PC3, data = Plasma_Proteins_Jan, colour = Subject) +
            scale_colour_manual(values = colorblind_Palette) 
pp3 <- qplot(x = PC1, y = PC2, data = Plasma_Proteins_Jan, colour = Subject) +
            scale_colour_manual(values = colorblind_Palette)
grid.arrange(pp1, pp2, pp3, ncol = 3,
             top = "Principal Components for Plasma Proteins(Visit 1)")
dev.off()

Plasma_Metabolites_Jan <- 
  read.csv(paste("chronobiome/Seth_PCA/PC_Scores/ScoresJan17/",
                 "PlasmaMetabolites_BothVisits_ScoresFor5PCs.csv"))
colnames(Plasma_Metabolites_Jan) <- c("Visit", "Subject", "Hour", "PC1",
                                      "PC2", "PC3", "PC4", "PC5")
postscript("Plasma_Metabs_PCA_Jan.eps")
pm1 <- qplot(x = PC2, y = PC3, data = Plasma_Metabolites_Jan, colour = Subject) +
            scale_colour_manual(values = colorblind_Palette)
pm2 <- qplot(x = PC1, y = PC3, data = Plasma_Metabolites_Jan, colour = Subject) +
            scale_colour_manual(values = colorblind_Palette) 
pm3 <- qplot(x = PC1, y = PC2, data = Plasma_Metabolites_Jan, colour = Subject) +
            scale_colour_manual(values = colorblind_Palette)
grid.arrange(pm1, pm2, pm3, ncol = 3, top = paste("Principal Components for ",
                                              "Plasma Metabolites(Visits", 
                                              " 1 and 2)")
             )
dev.off()


# Use Full.with.energy to save scatterplots of every combination of variables
# within the 48 hour window (406 plots)

Full.with.energy.reduced <- subset(Full.with.energy,
                                   select = -c(X1, Days, Times.y,
                                               TimeSubjectIndex, Times.x))

# Get a list of all 29-choose-2 combinations of two variables
combos <-  gtools::combinations(n = 29, r = 2,
                                v = colnames(Full.with.energy.reduced),
                                repeats.allowed = F)

# Save indices in the 'combos' list from 1 to 406
sequence <- seq(1, dim(combos)[1])

# Break up combos list into groups of 4 (and one group of 2)
# Save to "scatterplots.pdf" with 4 of the 406 plots are on
# each page 
chunks <- split(sequence, ceiling(seq_along(sequence) / 4))

pdf("scatterplots.pdf")
for(chunk in chunks){
  if(length(chunk) == 4){
    pair1 <- combos[chunk[1], ]
    pair2 <- combos[chunk[2], ]
    pair3 <- combos[chunk[3], ]
    pair4 <- combos[chunk[4], ]
    
    plot1 <- qplot(x = Full.with.energy.reduced[pair1[1]],
                   y = Full.with.energy.reduced[pair1[2]],
                   data = Full.with.energy.reduced, xlab = pair1[1], ylab = pair1[2])
    plot2 <- qplot(x = Full.with.energy.reduced[pair2[1]],
                   y = Full.with.energy.reduced[pair2[2]],
                   data = Full.with.energy.reduced, xlab = pair2[1], ylab = pair2[2])
    plot3 <- qplot(x = Full.with.energy.reduced[pair3[1]],
                   y = Full.with.energy.reduced[pair3[2]],
                   data = Full.with.energy.reduced, xlab = pair3[1], ylab = pair3[2])
    plot4 <- qplot(x = Full.with.energy.reduced[pair4[1]],
                   y = Full.with.energy.reduced[pair4[2]],
                   data = Full.with.energy.reduced, xlab = pair4[1], ylab = pair4[2])
    grid.arrange(plot1, plot2, plot3,plot4, ncol = 2, nrow = 2)

  } else {
    plot1 <- qplot(x = Full.with.energy.reduced[pair1[1]],
                   y = Full.with.energy.reduced[pair1[2]],
                   data = Full.with.energy.reduced, xlab = pair1[1], ylab = pair1[2])
    plot2 <- qplot(x = Full.with.energy.reduced[pair2[1]],
                   y = Full.with.energy.reduced[pair2[2]],
                   data = Full.with.energy.reduced, xlab = pair2[1], ylab = pair2[2])
    grid.arrange(plot1, plot2, ncol = 2, nrow = 2)
  }
}

dev.off()
