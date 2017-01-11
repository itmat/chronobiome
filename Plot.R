################################################
# Amy Campbell 2016
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


PlotR2s <- function(dataframe, title){
  # Generates geom_tile plot of the % variability explained (R2) for each variable by each other variable
  # :param: dataframe - Dataframe to plot (each variable in the dataframe will be compared to each other
  # :param: title - title to display on plot 
  #  dataframe, so these must be numeric)
  r2Matrix <- c()
  for (a in colnames(dataframe)) {
    for (b in colnames(dataframe)) {
      r2Matrix <-
        rbind(r2Matrix,
              GetR2(as.matrix(dataframe[a]),
                    as.matrix(dataframe[b])
              ))
    }
  }
  PlotMatrix <- data.frame(r2Matrix)
  colnames(PlotMatrix) <- c("factor.1", "factor.2", "variability.explained")

  PlotMatrix$variability.explained <- as.numeric(as.character(PlotMatrix$variability.explained))
  
  # Order variables for 4 month activity/communication dataset
  if(title=="Variability Explained in Activity and Communication Data"){
    f1 = factor(PlotMatrix$factor.1, levels(PlotMatrix$factor.1)[c(7, 8, 9 , 6, 18, 22, 21, 20, 4, 16, 15, 14, 1, 2, 3, 5, 17, 19, 10, 13, 12, 11)])
    f2 = factor(PlotMatrix$factor.2, levels(PlotMatrix$factor.2)[c(7, 8, 9 , 6, 18, 22, 21, 20, 4, 16, 15, 14, 1, 2, 3, 5, 17, 19, 10, 13, 12, 11)])           
  }
  # Order variables for 48 hour communication/activity/BP dataset
  else{
    f1 <- factor(PlotMatrix$factor.1, levels(PlotMatrix$factor.1)[c(8, 9, 10, 7, 21, 26, 25, 24, 5, 19, 18, 17, 1,
                                                                    2, 3, 6, 20, 22, 13, 16, 15, 14, 23, 4, 11, 27, 12)])
    f2 <- factor(PlotMatrix$factor.2, levels(PlotMatrix$factor.2)[c(8, 9, 10, 7, 21, 26, 25, 24, 5, 19, 18, 17, 1,
                                                                    2, 3, 6, 20, 22, 13, 16, 15, 14, 23, 4, 11, 27, 12)])
  }
  # generate geom_tile plot
  plot <-
    ggplot(data = PlotMatrix, aes(x = f1, y = f2, fill = variability.explained)) +
    geom_tile(color ="black") + scale_fill_gradient(low="white", high="steelblue4") +
    theme(axis.text.x = element_text(angle = 90, size =11), axis.text.y = element_text(size = 11)) +
    ggtitle(title) + xlab("Factor 1") + ylab("Factor 2")
  return(plot)
}

parsesubject <- function(row, half){
  # Function to be applied to all rows of 
  # the dataframe in question
  # :param: row - the entire row of the dataframe
  # :param: half - indicates whether the "subject" information
  # is in the first or second half of the string
  timesubjectindex <- row[1]
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
window1.begin <- as.numeric((chron("12/01/2014", format=c(dates="m/d/y")) - chron("2014-10-21", format=c(dates="y-m-d")) -1)*24 + 10)
window1.end <- as.numeric((chron("12/03/2014", format=c(dates="m/d/y")) - chron("2014-10-21", format=c(dates="y-m-d"))-1)*24 + 10)

window2.begin <- as.numeric((chron("12/15/2014", format=c(dates="m/d/y")) - chron("2014-10-21", format=c(dates="y-m-d")) -1)*24 + 10)
window2.end <- as.numeric((chron("12/17/2014", format=c(dates="m/d/y")) - chron("2014-10-21", format=c(dates="y-m-d"))-1)*24 + 10)

TimeList<- c(seq(window1.begin, window1.end), seq(window2.begin, window2.end))


###############################
# Heatmap of variance explained
###############################

transformed.communication.activity <- read.csv("transformed.communication.activity.csv")
transformed.communication.activity["X"] = NULL

# 48 Hours containing measurements of blood pressure, heart rate
BP_HR <- read.csv("heartrate.bp.csv")
Full <- dplyr::full_join(BP_HR, transformed.communication.activity, by="TimeSubjectIndex")
bp.HR.com.act <- na.omit(Full)

postscript("Variability.Act.Com.BP.eps")
PlotR2s(bp.HR.com.act[,5:dim(bp.HR.com.act)[2]], "Variability Explained in Activity, Communication, Biometric Data (48 Hours)")
dev.off()

# All four months (no blood pressure/heart rate data)
transformed.communication.activity["TimeSubjectIndex"] <- NULL
postscript("Variability.Activity.Communication.eps",  width = 480, height = 480)
PlotR2s(transformed.communication.activity, "Variability Explained in Activity and Communication Data")
dev.off()

################################
# Communication and activity PCA
################################

communication.activity <- (read.csv("communication.activity.csv"))
communication.activity['subject'] <- (apply(communication.activity, 1, parsesubject, 2))

communication.activity <- communication.activity[c()]
subj <- data.frame(communication.activity$subject)
length <- as.numeric(dim(communication.activity)[2])

Matrix_For_PCA <- communication.activity[,3:(length-1)]
Act_Com_PCA <- prcomp(Matrix_For_PCA, scale.=TRUE)
ActCom_PCA_Matrix <- data.frame(Act_Com_PCA$x)

ActCom_PCA_Matrix["Subject"] <- subj
write.table(Act_Com_PCA$rotation, "Act_Com_Loadings.csv", sep=",")

# Subsetting for 48 Hour Window
window1.begin <- as.numeric((chron("12/01/2014", format=c(dates="m/d/y")) - chron("2014-10-21", format=c(dates="y-m-d")) -1)*24 + 10)
window1.end <- as.numeric((chron("12/03/2014", format=c(dates="m/d/y")) - chron("2014-10-21", format=c(dates="y-m-d"))-1)*24 + 10)

window2.begin <- as.numeric((chron("12/15/2014", format=c(dates="m/d/y")) - chron("2014-10-21", format=c(dates="y-m-d")) -1)*24 + 10)
window2.end <- as.numeric((chron("12/17/2014", format=c(dates="m/d/y")) - chron("2014-10-21", format=c(dates="y-m-d"))-1)*24 + 10)

TimeList<- c(seq(window1.begin, window1.end), seq(window2.begin, window2.end))

TimeSubj <- data.frame(rownames(ActCom_PCA_Matrix))
ActCom_PCA_Matrix["TimesOnly"] <- apply(TimeSubj, 1, parsesubject, 1)

windowed_ActCom_PCA <- subset(ActCom_PCA_Matrix, TimesOnly %in% TimeList)

###########
# PCA Plots
###########
colorblind_Palette <- c("#000000", "#0072B2", "#56B4E9", "#F0E442", "#D55E00", "#CC79A7")

postscript("Act_Com_PCAPlotAll.eps", width=480, height=480)
p1 <- qplot(x=PC2, y=PC3, data=ActCom_PCA_Matrix, colour=Subject) + scale_colour_manual(values=colorblind_Palette)
p2 <- qplot(x=PC1, y=PC3, data=ActCom_PCA_Matrix, colour=Subject) + scale_colour_manual(values=colorblind_Palette)
p3 <- qplot(x=PC1, y=PC2, data=ActCom_PCA_Matrix, colour=Subject) + scale_colour_manual(values=colorblind_Palette)
grid.arrange(p1, p2, p3, ncol=3, top="Principal Components for Combined Activity and Communication (All Measurements)")
dev.off()

postscript("Act_Com_PCAPlot_48.eps")
p_1 <- qplot(x=PC2, y=PC3, data=windowed_ActCom_PCA, colour=Subject) + scale_colour_manual(values=colorblind_Palette) 
p_2 <- qplot(x=PC1, y=PC3, data=windowed_ActCom_PCA, colour=Subject) + scale_colour_manual(values=colorblind_Palette)
p_3 <- qplot(x=PC1, y=PC2, data=windowed_ActCom_PCA, colour=Subject) + scale_colour_manual(values=colorblind_Palette)
grid.arrange(p_1, p_2, p_3, ncol=3, top="Principal Components for Combined Activity and Communication (Visits 1 and 2)")
dev.off()

postscript("Facets_Act_Com_PC2PC3.eps")
fp_1 <- p_1 + facet_grid(. ~Subject)
grid.arrange(fp_1, ncol=1, top ="Principal Components for Combined Activity and Communication (Visits 1 and 2)")
dev.off()

postscript("Facets_Act_Com_PC1PC3.eps")
fp_2 <- p_2 + facet_grid(. ~Subject) 
grid.arrange(fp_2, ncol=1, top ="Principal Components for Combined Activity and Communication (Visits 1 and 2)")
dev.off()

postscript("Facets_Act_Com_PC1PC2.eps")
fp_3 <- p_3 + facet_grid(. ~Subject) 
grid.arrange(fp_3, ncol=1, top ="Principal Components for Combined Activity and Communication (Visits 1 and 2)")
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

