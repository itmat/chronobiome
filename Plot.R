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
  # Generates geom_tile plot of the % variance explained (R2) for each variable by each other variable
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
  # :param: Sequence - string indicating which order of ints to use to plot variables
  #                    on the variance explained output plot:  
  #                    sequence_actcom (just activity & communication variables)
  #                    sequence_actcomBP (activity & communication with blood pressure/HR data)
  #                    sequence_energy (for datasets including activity, communication,
  #                                     blood pressure, and energy data)
              
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
  PlotMatrix$variability.explained <- as.numeric(as.character(PlotMatrix$variability.explained))

  # Order variables using their alphabetically-determined level indices such that activraphy, actigraphy circadian, 
  # communication, communication circadian, and biometric variables are grouped together in the output chart
  sequence_actcomBP = c(8, 9, 10, 7, 21, 26, 25, 24, 5, 19, 18, 17, 1,
                        2, 3, 6, 20, 22, 13, 16, 15, 14, 23, 4, 11, 27, 12)
  sequence_actcom =  c(7, 8, 9 , 6, 18, 22, 21, 20, 4, 16, 15, 14, 1,
                       2, 3, 5, 17, 19, 10, 13, 12, 11)
  sequence_energy = c(8, 9, 10, 7, 23, 28, 27, 26, 5, 21, 20, 18, 1,
                      2, 3, 6, 22, 24, 13, 16, 15, 14, 25, 4, 11, 29, 12, 17, 19)
  
  if(Sequence == "ActCom"){
    f1 = factor(PlotMatrix$factor.1, levels(PlotMatrix$factor.1)[sequence_actcom])
    f2 = factor(PlotMatrix$factor.2, levels(PlotMatrix$factor.2)[sequence_actcom])           
  }
  # Order variables for 48 hour communication/activity/BP dataset
  else if(Sequence=="BP"){
    f1 <- factor(PlotMatrix$factor.1, levels(PlotMatrix$factor.1)[sequence_actcomBP])
    f2 <- factor(PlotMatrix$factor.2, levels(PlotMatrix$factor.2)[sequence_actcomBP])
  }
  else{
    f1 <- factor(PlotMatrix$factor.1, levels(PlotMatrix$factor.1)[sequence_energy])
    f2 <- factor(PlotMatrix$factor.2, levels(PlotMatrix$factor.2)[sequence_energy])
  }
  # generate geom_tile plot
  plot <-
    ggplot(data = PlotMatrix, aes(x = f1, y = f2, fill = variability.explained)) +
    geom_tile(color ="black") + scale_fill_gradient(low="white", high="steelblue4") +
    theme(axis.text.x = element_text(angle = 90, size =10), axis.text.y = element_text(size = 10)) +
    ggtitle(title) + xlab("Factor 1") + ylab("Factor 2")
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
  
  if(toString(half)=="2"){
    return(unlist(strlist)[2])
  }
  else{
    return(as.numeric(unlist(strlist)[1]))
  }
}

############################################
# List of times to include in 48 hour window
############################################
window1.begin <- as.numeric((chron("12/01/2014", format=c(dates="m/d/y")) - chron("2014-10-21",
                                                          format=c(dates="y-m-d")) -1)*24 + 10)
window1.end <- as.numeric((chron("12/03/2014", format=c(dates="m/d/y")) - chron("2014-10-21", 
                                                          format=c(dates="y-m-d"))-1)*24 + 10)

window2.begin <- as.numeric((chron("12/15/2014", format=c(dates="m/d/y")) - chron("2014-10-21", 
                                                          format=c(dates="y-m-d")) -1)*24 + 10)

window2.end <- as.numeric((chron("12/17/2014", format=c(dates="m/d/y")) - chron("2014-10-21", 
                                                          format=c(dates="y-m-d"))-1)*24 + 10)

TimeList<- c(seq(window1.begin, window1.end), seq(window2.begin, window2.end))
TimeList_window1 <- c(seq(window1.begin, window1.end))
TimeList_window2 <- c(seq(window2.begin, window2.end))

###############################
# Heatmap of variance explained
###############################

transformed.communication.activity <- readr::read_csv("transformed.communication.activity.csv")
transformed.communication.activity["X1"] = NULL

# Two 48 Hour visits containing measurements of blood pressure, heart rate
BP_HR <- readr::read_csv("heartrate.bp.csv")
Full <- dplyr::full_join(BP_HR, transformed.communication.activity, by="TimeSubjectIndex")
bp.HR.com.act <- na.omit(Full)
postscript("Variability.Act.Com.BP.eps")
PlotR2s(bp.HR.com.act[,5:dim(bp.HR.com.act)[2]], 
        "Variance Explained in Activity, Communication, Biometric Data (Visits 1 and 2)", "BP")
dev.off()

# 48 hour visit 1
Visit1_bp.HR.com.act <- subset(bp.HR.com.act, Times %in% TimeList_window1)

postscript("Varability.Act.Com.BP.Visit1")
PlotR2s(Visit1_bp.HR.com.act[, 5:dim(Visit1_bp.HR.com.act)[2]],
        "Variance Explained in Activity, Communication, Biometric Data(Visit 1)", "BP")
dev.off()

# 48 hour visit 2
Visit2_bp.HR.com.act <- subset(bp.HR.com.act, Times %in% TimeList_window2)
postscript("Varability.Act.Com.BP.Visit2")
PlotR2s(Visit2_bp.HR.com.act[, 5:dim(Visit2_bp.HR.com.act)[2]],
        "Variance Explained in Activity, Communication, Biometric Data (Visit 2)", "BP")
dev.off()


# All four months (no blood pressure/heart rate data)
transformed.communication.activity["TimeSubjectIndex"] <- NULL
postscript("Variability.Activity.Communication.eps",  width = 480, height = 480)
PlotR2s(transformed.communication.activity,
        "Variance Explained in Activity and Communication Data", "ActCom")
dev.off()

# Both Visits by Subject 
bp.HR.com.act["Subject"] <- (apply(bp.HR.com.act, 1, parsesubject, 2))

# Subject HCR001 
HCR001_set <- subset(bp.HR.com.act, Subject=="HCR001")
HCR001_set["Subject"] <- NULL
postscript("HCR001_Variance.eps",  width = 480, height = 480)
PlotR2s(HCR001_set[, 5:dim(HCR001_set)[2]], "Variance Explained in HCR001", "BP")
dev.off()

# HCR003
HCR003_set <- subset(bp.HR.com.act, Subject=="HCR003")
HCR003_set["Subject"] <- NULL
postscript("HCR003_Variance.eps",  width = 480, height = 480)
PlotR2s(HCR003_set[, 5:dim(HCR003_set)[2]], "Variance Explained in HCR003", "BP")
dev.off()


HCR004_set <- subset(bp.HR.com.act, Subject=="HCR004")
HCR004_set["Subject"] <- NULL
postscript("HCR004_Variance.eps",  width = 480, height = 480)
PlotR2s(HCR004_set[, 5:dim(HCR004_set)[2]], "Variance Explained in HCR004", "BP")
dev.off()

# HCR006
HCR006_set <- subset(bp.HR.com.act, Subject=="HCR006")
HCR006_set["Subject"] <- NULL
postscript("HCR006_Variance.eps",  width = 480, height = 480)
PlotR2s(HCR006_set[, 5:dim(HCR006_set)[2]], "Variance Explained in HCR006", "BP")
dev.off()

# HCR008
HCR008_set <- subset(bp.HR.com.act, Subject=="HCR008")
HCR008_set["Subject"] <- NULL
postscript("HCR008_Variance.eps",  width = 480, height = 480)
PlotR2s(HCR008_set[, 5:dim(HCR008_set)[2]], "Variance Explained in HCR008", "BP")
dev.off()

# HCR009
HCR009_set <- subset(bp.HR.com.act, Subject=="HCR009")
HCR009_set["Subject"] <- NULL
postscript("HCR009_Variance.eps",  width = 480, height = 480)
PlotR2s(HCR009_set[, 5:dim(HCR009_set)[2]], "Variance Explained in HCR009", "BP")
dev.off()

energy <- read.csv("energy.csv")
energy["X"] <- NULL
energy["Times"] <- apply(energy, 1, parsesubject, 1)

Full.with.energy <- dplyr::full_join(bp.HR.com.act, energy, by="TimeSubjectIndex")
Full.with.energy <- na.omit(Full.with.energy)
Full.with.energy["Subject"] <- NULL
postscript("Variance_Explained_WithEnergy.eps")
PlotR2s(Full.with.energy[, 6:dim(Full.with.energy)[2]-1], "Variance Explained in Activity, Communication, Blood Pressure, and Energy Variables", "somethingelse")
dev.off()
################################
# Communication and activity PCA
################################

communication.activity <- (read.csv("communication.activity.csv"))
communication.activity['subject'] <- (apply(communication.activity, 1, parsesubject, 2))

TimeSubj <- data.frame(communication.activity$TimeSubjectIndex)
subj <- data.frame(communication.activity$subject)
length <- as.numeric(dim(communication.activity)[2])
Matrix_For_PCA <- communication.activity[,3:(length-1)]
Act_Com_PCA <- prcomp(Matrix_For_PCA, scale.=TRUE)
ActCom_PCA_Matrix <- data.frame(Act_Com_PCA$x)

ActCom_PCA_Matrix["Subject"] <- subj
write.table(Act_Com_PCA$rotation, "Act_Com_Loadings.csv", sep=",")


ActCom_PCA_Matrix["TimeSubjectIndex"] <- TimeSubj
  
ActCom_PCA_Matrix["TimesOnly"] <- apply(ActCom_PCA_Matrix, 1, parsesubject, 1)

windowed_ActCom_PCA <- subset(ActCom_PCA_Matrix, TimesOnly %in% TimeList)

###########
# PCA Plots
###########
colorblind_Palette <- c("#000000", "#0072B2", "#56B4E9", "#F0E442", "#D55E00", "#CC79A7")

postscript("Act_Com_PCAPlotAll.eps", width=480, height=480)
p1 <- qplot(x=PC2, y=PC3, data=ActCom_PCA_Matrix, colour=Subject) + 
            scale_colour_manual(values=colorblind_Palette)
p2 <- qplot(x=PC1, y=PC3, data=ActCom_PCA_Matrix, colour=Subject) + 
            scale_colour_manual(values=colorblind_Palette)
p3 <- qplot(x=PC1, y=PC2, data=ActCom_PCA_Matrix, colour=Subject) + 
            scale_colour_manual(values=colorblind_Palette)
grid.arrange(p1, p2, p3, ncol=3, 
             top="Principal Components for Combined Activity and Communication (All Measurements)")
dev.off()

postscript("Act_Com_PCAPlot_48.eps")
p_1 <- qplot(x=PC2, y=PC3, data=windowed_ActCom_PCA, colour=Subject) +
             scale_colour_manual(values=colorblind_Palette) 
p_2 <- qplot(x=PC1, y=PC3, data=windowed_ActCom_PCA, colour=Subject) +
             scale_colour_manual(values=colorblind_Palette)
p_3 <- qplot(x=PC1, y=PC2, data=windowed_ActCom_PCA, colour=Subject) +
             scale_colour_manual(values=colorblind_Palette)
grid.arrange(p_1, p_2, p_3, ncol=3,
             top="Principal Components for Combined Activity and Communication (Visits 1 and 2)")
dev.off()

postscript("Facets_Act_Com_PC2PC3.eps")
fp_1 <- p_1 + facet_grid(. ~Subject)
grid.arrange(fp_1, ncol=1,
             top ="Principal Components for Combined Activity and Communication (Visits 1 and 2)")
dev.off()

postscript("Facets_Act_Com_PC1PC3.eps")
fp_2 <- p_2 + facet_grid(. ~Subject) 
grid.arrange(fp_2, ncol=1, 
             top ="Principal Components for Combined Activity and Communication (Visits 1 and 2)")
dev.off()

postscript("Facets_Act_Com_PC1PC2.eps")
fp_3 <- p_3 + facet_grid(. ~Subject) 
grid.arrange(fp_3, ncol=1, 
             top ="Principal Components for Combined Activity and Communication (Visits 1 and 2)")
dev.off()


# Load and plot other PCAs (Generated by S. Rhodes)


Saliva_Metabolites <- read.csv("chronobiome/Seth_PCA/PC_Scores/SalivaMetabolites_BothVisits_ScoresFor10PCs.csv")
colnames(Saliva_Metabolites) <- c("Visit", "Subject", "Hour", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
postscript("Saliva_Metabs_PCA.eps")
sm1 <- qplot(x=PC2, y=PC3, data=Saliva_Metabolites, colour=Subject) + scale_colour_manual(values=colorblind_Palette)
sm2 <- qplot(x=PC1, y=PC3, data=Saliva_Metabolites, colour=Subject) + scale_colour_manual(values=colorblind_Palette)
sm3 <- qplot(x=PC1, y=PC2, data=Saliva_Metabolites, colour=Subject) + scale_colour_manual(values=colorblind_Palette)
grid.arrange(sm1, sm2, sm3, ncol=3, top="Principal Components for Saliva Metabolites(Visits 1 and 2)")
dev.off()

Saliva_Microbes <- read.csv("chronobiome/Seth_PCA/PC_Scores/SalivaMicrobes_BothVisits_ScoresFor5PCs.csv")
colnames(Saliva_Microbes) <-(c("X", "Subject", "Hour", "PC1", "PC2", "PC3", "PC4", "PC5"))
postscript("Saliva_Microbes_PCA.eps")
smic1 <- qplot(x=PC2, y=PC3, data=Saliva_Microbes, colour=Subject) + scale_colour_manual(values=colorblind_Palette) 
smic2 <- qplot(x=PC1, y=PC3, data=Saliva_Microbes, colour=Subject) + scale_colour_manual(values=colorblind_Palette) 
smic3 <- qplot(x=PC1, y=PC2, data=Saliva_Microbes, colour=Subject) + scale_colour_manual(values=colorblind_Palette)
grid.arrange(smic1, smic2, smic3, ncol=3, top="Principal Components for Saliva Microbes(Visits 1 and 2)")
dev.off()

Plasma_Proteins <- read.csv("chronobiome/Seth_PCA/PC_Scores/PlasmaProteins_1stVisits_HCR009-12hrRemove_ScoresFor10PCs.csv")
Plasma_Proteins <- na.omit(Plasma_Proteins)
colnames(Plasma_Proteins) <- c("X", "Subject", "Hour", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7","PC8","PC9","PC10")
postscript("Plasma_Proteins_PCA.eps")
pp1 <- qplot(x=PC2, y=PC3, data=Plasma_Proteins, colour=Subject) + scale_colour_manual(values=colorblind_Palette) 
pp2 <- qplot(x=PC1, y=PC3, data=Plasma_Proteins, colour=Subject) + scale_colour_manual(values=colorblind_Palette) 
pp3 <- qplot(x=PC1, y=PC2, data=Plasma_Proteins, colour=Subject) + scale_colour_manual(values=colorblind_Palette)
grid.arrange(pp1, pp2, pp3, ncol=3, top="Principal Components for Plasma Proteins(Visit 1)")
dev.off()

Plasma_Metabolites <- read.csv("chronobiome/Seth_PCA/PC_Scores/PlasmaMetabolites_BothVisits_ScoresFor10PCs.csv")
colnames(Plasma_Metabolites) <- c("Visit", "Subject", "Hour", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
postscript("Plasma_Metabs_PCA.eps")
pm1 <- qplot(x=PC2, y=PC3, data=Plasma_Metabolites, colour=Subject) + scale_colour_manual(values=colorblind_Palette)
pm2 <- qplot(x=PC1, y=PC3, data=Plasma_Metabolites, colour=Subject) + scale_colour_manual(values=colorblind_Palette) 
pm3 <- qplot(x=PC1, y=PC2, data=Plasma_Metabolites, colour=Subject) + scale_colour_manual(values=colorblind_Palette)
grid.arrange(pm1, pm2, pm3, ncol=3, top="Principal Components for Plasma Metabolites(Visits 1 and 2)")
dev.off()

Macronutrients <- read.csv("/home/amycampbell/Downloads/MacroNutrientData_SmallerChunks_BothVisits.csv")
postscript("Macronutrients.eps")
n1 <- qplot(x=Protein, y=TotalFat, data=Macronutrients, colour=Subject) + scale_colour_manual(values=colorblind_Palette)
n2 <- qplot(x=Protein, y=Carbohydrate, data=Macronutrients, colour=Subject) + scale_colour_manual(values=colorblind_Palette)
n3 <- qplot(x=Carbohydrate, y=TotalFat, data=Macronutrients, colour=Subject) + scale_colour_manual(values=colorblind_Palette)
grid.arrange(n1, n2, n3, ncol=3, top="Macronutrients(Visits 1 and 2)")
dev.off()

