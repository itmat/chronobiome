################################################
# Amy Campbell 2016-2017
#
# Plots variability explained by variables in activity, communication,
# energy, and blood pressure/heart rate data 
# 
# Cluster subjects by lm-fit p-values and draw dendrograms
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

`%>%` <- magrittr::`%>%`

##################
# Define functions
##################


GetR2 <- function(var1, var2, return_pvalue=FALSE) {
  # Generates a simple linear regression model of var1~var2, and returns
  # R-squared for that model to represent proportion of V1 explained by V2
  # :param: var1 - Variable to 'be explained' by var2
  # :param: var2 - Variable to 'explain' var1
  # :param: return_pvalue - boolean indicating whether or not function returns
  #                         p-value for the fit of the estimated linear model
  #                         to the data, instead of the R-squared value:
  #                         FALSE (default) - return R-sqaured value.
  #                         TRUE - return p-value from fit.
  Var1name <- colnames(var1)[1]
  Var2name <-  colnames(var2)[1]
  model <- lm(as.numeric(var1) ~ as.numeric(var2))
  
  ret_value = summary(model)$r.squared
  
  if(return_pvalue) {
      #Calculate the p-value using F-statistic attributes from linear model
      fstats <- summary(model)$fstatistic
      
      ret_value <- NA
      
      #If an f-test was performed successfully
      if(!is.null(fstats)) {
          pvalue <- pf(fstats[1],fstats[2],fstats[3],lower.tail=F)
          attributes(pvalue) <- NULL
          ret_value <- pvalue
      }
      
  }
  
  return(c(Var1name, Var2name, ret_value))
}


PlotR2s <- function(chrono_df, title, Sequence, triangle_heatmap=FALSE, add_labels=FALSE, color_by_pvalue=FALSE) {
  # Generates geom_tile plot of the % variance explained (R2) for each variable
  # by each other variable.
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
  # :param: triangle_heatmap - boolean indicating whether to plot just the
  #                            top of the heatmap, above the diagonal. The
  #                            heatmap is mirrored across the diagonal, so
  #                            plotting the full heatmap technically contains
  #                            redundant information.
  # :param: add_labels - boolean indicating whether to add text labels to each
  #                      box in the heatmap that list the R-squared values.
  #                      The listed values are rounded to two decimal places.
  # :param: color_by_pvalue - boolean indicating whether to color the heatmap
  #                           using p-values from the fit of the lm model,
  #                           instead of the R-squared values:
  #                           FALSE (default) - Color and label plot using
  #                                             R-squared values.
  #                           TRUE - Color and label plot using p-values.
  #                           Note, this parameter also causes the graph to
  #                           display p-values instead of R-squared values
  #                           when add_labels is set to TRUE.
    
  chrono_df["Subject"] <- NULL
  chrono_df["TimeSubjectIndex"] <- NULL            
  r2Matrix <- c()
  for (a in colnames(chrono_df)) {
    for (b in colnames(chrono_df)) {
      r2Matrix <-
        rbind(r2Matrix,
              GetR2(as.matrix(chrono_df[a]),
                    as.matrix(chrono_df[b]),
                    color_by_pvalue
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
  
  #Update factors with new ordering
  plot.matrix$factor.1 = f1
  plot.matrix$factor.2 = f2
  
  if(triangle_heatmap) {
      #Remove data from lower triangle
      plot.matrix <- dcast(plot.matrix, formula = factor.1 ~ factor.2, value.var = "variability.explained")
      #Exclude factor.1 labels when determining triangle plot
      factor1 = plot.matrix$factor.1
      plot.matrix$factor.1 <- NULL
      plot.matrix[lower.tri(plot.matrix)] <- NA
      plot.matrix$factor.1 = factor1
      plot.matrix <- melt(plot.matrix,  variable.name = "factor.2",
                          value.name = "variability.explained",
                          id.vars = c("factor.1"))
      plot.matrix = na.omit(plot.matrix)
  }
  
  # generate geom_tile plot
  plot <-
    ggplot(data = plot.matrix,
           aes(x = factor.1, y = factor.2, fill = variability.explained)) +
           geom_tile(color = "black")
  
  # Adjust color scale for displaying either R-squared values or p-values.
  if(color_by_pvalue) { #Color by p-values
      plot <- plot + scale_fill_gradientn(name="p-value",
                                          colors = c("steelblue4","white","white"),
                                          values = c(0,0.05,1))
  } else { #Color by R-squared
      plot <- plot + scale_fill_gradientn(name="variability\nexplained",
                                          colors = c("white", "steelblue4"),
                                          values = c(0,1))
  }
  
  plot <-
      plot +
      theme(axis.text.x = element_text(angle = 90, size = 8,
                                       hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 8),
            panel.background = element_blank()) +
      ggtitle(title) +
      xlab("Factor 1") + ylab("Factor 2")
  
  if(add_labels) {
      plot <- plot + geom_text(aes(label=sprintf("%0.2f", round(variability.explained, digits = 2))),
                               size = 1.75)
  }

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


FormatPvalueTableForClustering <- function(pvalue_table) {
    # Given a p-value table, remove self comparisons
    # (where factor.1 == factor.2) and combine
    # reciprocal comparisons (e.g. factor.1=log.Steps
    # and factor.2=log.Axis1 vs. factor.1=log.Axis1
    # and factor.2=log.Steps).
    # :param: pvalue_table - dataframe with the following
    #                        columns: factor.1, factor.2,
    #                        and variability.explained)
    
    formatted.pvalue_table <-
        pvalue_table %>% 
        dplyr::rename(p_value = variability.explained) %>% 
        #Remove diagonal, since fit is perfect
        dplyr::filter(as.character(factor.1) != as.character(factor.2)) %>% 
        #Combine factor.1 and factor.2, sorting each individual pair in alphabetical order.
        #This makes it easy to remove reciprocal comparisons (e.g. log.Steps_log.Axis1 and
        #log.Axis1_log.Steps).
        dplyr::mutate(comparison.label = ifelse(as.character(factor.1) < as.character(factor.2),
                                                paste0(as.character(factor.1), "_", as.character(factor.2)),
                                                paste0(as.character(factor.2), "_", as.character(factor.1))),
                      comparison.order = ifelse(as.character(factor.1) < as.character(factor.2),
                                                1,2)) %>% 
        dplyr::select(comparison.label, comparison.order, p_value) %>%
        dplyr::arrange(comparison.label, comparison.order) %>% 
        #Due to floating-point errors, there can be slight differences between the p-values
        #generated by the reciprocal comparisons. While these differences are negligible,
        #they are enough to prevent filtering using the distinct() and unique() functions.
        #To overcome this, the code below separates the p-values form each reciprocal
        #comparison into two separate columns, and only retains the smaller of the two.
        dcast(formula = comparison.label ~ comparison.order, value.var = "p_value") %>% 
        dplyr::mutate(p_value = pmin(`1`, `2`)) %>% 
        dplyr::select(comparison.label, p_value)
    
    return(formatted.pvalue_table)
}


PlotDendro <- function(pvalue_table, title) {
    # Plots a dendrogram, given a formatted table of
    # LM-fit p-values. This version of the function
    # uses base R graphics.
    # :param: pvalue_table - Dataframe with the following columns:
    #                        Subject, comparison.label, p_value
    # :param: title - title to display on plot
    
    data_for_clustering <- 
        pvalue_table %>% 
        dcast(formula = Subject ~ comparison.label, value.var = "p_value") %>% 
        tibble::column_to_rownames(var="Subject") %>% 
        as.matrix()
    
    plot(hclust(dist(data_for_clustering)),
         main=title,
         sub="",
         xlab="",
         ylab="Euclidean distance")
}


PlotDendro.w_ggdendro <- function(pvalue_table, title) {
    # Generates a dendrogram as a ggplot object, given a
    # formatted table of LM-fit p-values. This function
    # is appropriate if the user needs to pass the plot
    # as an object, or wishes to leverage ggplot2 function
    # and features to further modify the plot.
    # :param: pvalue_table - Dataframe with the following columns:
    #                        Subject, comparison.label, p_value
    # :param: title - title to display on plot
    
    data_for_clustering <- 
        pvalue_table %>% 
        dcast(formula = Subject ~ comparison.label, value.var = "p_value") %>% 
        tibble::column_to_rownames(var="Subject") %>% 
        as.matrix()
    
    data_for_ggdendro <-
        ggdendro::dendro_data(as.dendrogram(hclust(dist(data_for_clustering)),
                                            type="rectangle"))
    
    plot <-
        ggplot(ggdendro::segment(data_for_ggdendro)) +
        geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
        scale_x_continuous(breaks = seq_along(data_for_ggdendro$labels$label), 
                           labels = data_for_ggdendro$labels$label) +
        ylab("Euclidean distance") +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              axis.title.x = element_blank(),
              axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1),
              axis.title.y = element_text(angle = 90, hjust = 0.5),
              panel.background = element_blank(),
              axis.ticks = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    
    return(plot)
}



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
              "Communication, Biometric Data (Visits 1 and 2)"), "BP",
        triangle_heatmap = FALSE, add_labels = FALSE)
dev.off()


#Define spans for visits 1 & 2 in terms of the standardized time units
TimeList_window1 <- c(seq(922, 1232))
TimeList_window2 <- c(seq(1328, 1545))

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
postscript("HCR004_Variance.4months_withEnergy.eps",
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
postscript("HCR001_Variance_visit1.with_labels.eps",  width = 480, height = 480)
PlotR2s(HCR001.set1[, 6:dim(HCR001.set1)[2] - 2],
        "Variance Explained in HCR001 (visit 1)", "energy",
        triangle_heatmap = TRUE, add_labels = TRUE)
dev.off()


# HCR003
HCR003.set1 <- subset(subset(Full.with.energy, Subject == "HCR003"), (Days <= 43))
postscript("HCR003_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR003.set1[, 6:dim(HCR003.set1)[2] - 2],
        "Variance Explained in HCR003 (visit 1)", "energy")
dev.off()
postscript("HCR003_Variance_visit1.with_labels.eps",  width = 480, height = 480)
PlotR2s(HCR003.set1[, 6:dim(HCR003.set1)[2] - 2],
        "Variance Explained in HCR003 (visit 1)", "energy",
        triangle_heatmap = TRUE, add_labels = TRUE)
dev.off()

# HCR004
HCR004.set1 <- subset(subset(Full.with.energy, Subject == "HCR004"), (Days<=45))
postscript("HCR004_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR004.set1[, 6:dim(HCR004.set1)[2] - 2],
        "Variance Explained in HCR004 (visit 1)", "energy")
dev.off()
postscript("HCR004_Variance_visit1.with_labels.eps",  width = 480, height = 480)
PlotR2s(HCR004.set1[, 6:dim(HCR004.set1)[2] - 2],
        "Variance Explained in HCR004 (visit 1)", "energy",
        triangle_heatmap = TRUE, add_labels = TRUE)
dev.off()

# HCR006
HCR006.set1 <- subset(subset(Full.with.energy, Subject == "HCR006"), (Days<=44))
postscript("HCR006_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR006.set1[, 6:dim(HCR006.set1)[2] - 2],
        "Variance Explained in HCR006 (visit 1)", "energy")
dev.off()
postscript("HCR006_Variance_visit1.with_labels.eps",  width = 480, height = 480)
PlotR2s(HCR006.set1[, 6:dim(HCR006.set1)[2] - 2],
        "Variance Explained in HCR006 (visit 1)", "energy",
        triangle_heatmap = TRUE, add_labels = TRUE)
dev.off()

# HCR008
HCR008.set1 <- subset(subset(Full.with.energy, Subject == "HCR008"), Days<=45)
postscript("HCR008_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR008.set1[, 6:dim(HCR008.set1)[2] - 2],
        "Variance Explained in HCR008 (visit 1)", "energy")
dev.off()
postscript("HCR008_Variance_visit1.with_labels.eps",  width = 480, height = 480)
PlotR2s(HCR008.set1[, 6:dim(HCR008.set1)[2] - 2],
        "Variance Explained in HCR008 (visit 1)", "energy",
        triangle_heatmap = TRUE, add_labels = TRUE)
dev.off()

# HCR009
HCR009.set1 <- subset(subset(Full.with.energy, Subject == "HCR009"), Days <= 51)
postscript("HCR009_Variance_visit1.eps",  width = 480, height = 480)
PlotR2s(HCR009.set1[, 6:dim(HCR009.set1)[2] - 2],
        "Variance Explained in HCR009 (visit 1)", "energy")
dev.off()
postscript("HCR009_Variance_visit1.with_labels.eps",  width = 480, height = 480)
PlotR2s(HCR009.set1[, 6:dim(HCR009.set1)[2] - 2],
        "Variance Explained in HCR009 (visit 1)", "energy",
        triangle_heatmap = TRUE, add_labels = TRUE)
dev.off()


# VISIT 2
# Subject HCR001 
HCR001.set2 <- subset(subset(Full.with.energy, Subject == "HCR001"), Days>=55)
postscript("HCR001_Variance_visit2.eps",  width = 480, height = 480)
PlotR2s(HCR001.set2[, 6:dim(HCR001.set2)[2] - 2],
        "Variance Explained in HCR001 (visit 2)", "energy")
dev.off()
postscript("HCR001_Variance_visit2.with_labels.eps",  width = 480, height = 480)
PlotR2s(HCR001.set2[, 6:dim(HCR001.set2)[2] - 2],
        "Variance Explained in HCR001 (visit 2)", "energy",
        triangle_heatmap = TRUE, add_labels = TRUE)
dev.off()

# HCR003
HCR003.set2 <- subset(subset(Full.with.energy, Subject == "HCR003"), Days >= 55)
postscript("HCR003_Variance_visit2.eps",  width = 480, height = 480)
PlotR2s(HCR003.set2[, 6:dim(HCR003.set2)[2] - 2],
        "Variance Explained in HCR003 (visit 2)", "energy")
dev.off()
postscript("HCR003_Variance_visit2.with_labels.eps",  width = 480, height = 480)
PlotR2s(HCR003.set2[, 6:dim(HCR003.set2)[2] - 2],
        "Variance Explained in HCR003 (visit 2)", "energy",
        triangle_heatmap = TRUE, add_labels = TRUE)
dev.off()

# HCR004
HCR004.set2 <- subset(subset(Full.with.energy, Subject == "HCR004"), Days >= 55)
postscript("HCR004_Variance_visit2.eps",  width = 480, height = 480)
PlotR2s(HCR004.set2[, 6:dim(HCR004.set2)[2] - 2],
        "Variance Explained in HCR004 (visit 2)", "energy")
dev.off()
postscript("HCR004_Variance_visit2.with_labels.eps",  width = 480, height = 480)
PlotR2s(HCR004.set2[, 6:dim(HCR004.set2)[2] - 2],
        "Variance Explained in HCR004 (visit 2)", "energy",
        triangle_heatmap = TRUE, add_labels = TRUE)
dev.off()

# HCR006
HCR006.set2 <- subset(subset(Full.with.energy, Subject == "HCR006"), Days >= 57)
postscript("HCR006_Variance_visit2.eps",  width = 480, height = 480)
PlotR2s(HCR006.set2[, 6:dim(HCR006.set2)[2] - 2],
        "Variance Explained in HCR006 (visit 2)", "energy")
dev.off()
postscript("HCR006_Variance_visit2.with_labels.eps",  width = 480, height = 480)
PlotR2s(HCR006.set2[, 6:dim(HCR006.set2)[2] - 2],
        "Variance Explained in HCR006 (visit 2)", "energy",
        triangle_heatmap = TRUE, add_labels = TRUE)
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
postscript("HCR009_Variance_visit2.with_labels.eps",  width = 480, height = 480)
PlotR2s(HCR009.set2[, 6:dim(HCR009.set2)[2] - 2],
        "Variance Explained in HCR009 (visit 2)", "energy",
        triangle_heatmap = TRUE, add_labels = TRUE)
dev.off()



# Heatmaps of LM-fit p-values & cluster subjects by LM-fit p-values
###################################################################

# Two 48 Hour visits containing measurements of blood pressure, heart rate
postscript("LM_pvalues.Act.Com.BP.eps")
PlotR2s(bp.HR.com.act[,5:dim(bp.HR.com.act)[2]], 
        paste0("Model fit p-values for Activity,", 
               "Communication, Biometric Data (Visits 1 and 2)"), "BP",
        triangle_heatmap = FALSE, add_labels = FALSE, color_by_pvalue = TRUE)
dev.off()

# 48 hour visit 1
postscript("LM_pvalues.Act.Com.BP.Visit1.Feb.eps")
PlotR2s(Visit1_bp.HR.com.act[, 5:dim(Visit1_bp.HR.com.act)[2]],
        "Model fit p-values for Activity, Communication, Biometric Data(Visit 1)",
        "BP", color_by_pvalue = TRUE)
dev.off()

# 48 hour visit 2
postscript("LM_pvalues.Act.Com.BP.Visit2.Feb.eps")
PlotR2s(Visit2_bp.HR.com.act[, 5:dim(Visit2_bp.HR.com.act)[2]],
        "Model fit p-values for Activity, Communication, Biometric Data (Visit 2)",
        "BP", color_by_pvalue = TRUE)
dev.off()

# Two 48 Hour visits containing measurements of blood pressure, heart rate, and energy expenditure
postscript("LM_pvalues_WithEnergy.eps")
PlotR2s(Full.with.energy[, 6:dim(Full.with.energy)[2] - 2],
        paste0("Model fit p-values for Activity, Communication, Blood Pressure,",
               " and Energy Variables"),
        "energy", color_by_pvalue = TRUE)
dev.off()

# All four months (no blood pressure/heart rate data)
postscript("LM_pvalues.Activity.Communication.eps",  width = 480, height = 480)
PlotR2s(transformed.communication.activity,
        "Model fit p-values for Activity and Communication Data", "ActCom",
        color_by_pvalue = TRUE)
dev.off()

# All four months (with energy)
postscript("LM_pvalues.Activity.Communication_WithEnergy.eps",
           width = 480, height = 480)
PlotR2s(Energy.4months,
        "Model fit p-values for Activity and Communication Data (4 months)",
        "4monthenergy", color_by_pvalue = TRUE)
dev.off()

#A note about clustering. The ggplot2 object returned by the PlotR2s function contains
#a table of the original data table used for plotting. The code below uses this fact
#to extract the table of LM-fit p-values for each subject while generating the heatmaps.


# All four months by subject (no blood pressure/heart rate data)
all_subjects.pvalue_table.4months = data.frame()

for(subject in c("HCR001","HCR003","HCR004","HCR006","HCR008","HCR009")) {
    subject.4months <- subset(transformed.communication.activity, Subject == subject)
    # postscript(paste0(subject, "_LM_pvalues.4months.eps"),  width = 480, height = 480)
    figure_plot = 
        PlotR2s(subject.4months, paste0("Model fit p-values for ", subject, " (all 4 months)"),
                "ActCom", color_by_pvalue = TRUE)
    # print(figure_plot)
    # dev.off()
    
    #Extract table of p-values for current subject
    subject.pvalue_table = figure_plot$data
    
    #Format and add to existing table of p-values for all subjects
    all_subjects.pvalue_table.4months <-
        FormatPvalueTableForClustering(subject.pvalue_table) %>% 
        dplyr::mutate(Subject = subject) %>% 
        dplyr::select(Subject, comparison.label, p_value) %>%
        dplyr::bind_rows(all_subjects.pvalue_table.4months)
}

#Cluster data by subjects and plot dendrogram
postscript("Clustering_dendro.LM_pvalues.Activity_Communication.4months.eps",  width = 480, height = 480)
PlotDendro(all_subjects.pvalue_table.4months,
           "Cluster subjects by LM-fit p-value using 4 months of data")
dev.off()


# All four months by subject (with energy)
all_subjects.pvalue_table.4months_w_E = data.frame()

for(subject in c("HCR001","HCR003","HCR004","HCR006","HCR008","HCR009")) {
    subject.4months <- subset(Energy.4months, Subject == subject)
    # postscript(paste0(subject, "_LM_pvalues.4months_withEnergy.eps"),  width = 480, height = 480)
    figure_plot = 
        PlotR2s(subject.4months, paste0("Model fit p-values for ", subject, " (all 4 months)"),
                "4monthenergy", color_by_pvalue = TRUE)
    # print(figure_plot)
    # dev.off()
    
    #Extract table of p-values for current subject
    subject.pvalue_table = figure_plot$data
    
    #Format and add to existing table of p-values for all subjects
    all_subjects.pvalue_table.4months_w_E <-
        FormatPvalueTableForClustering(subject.pvalue_table) %>% 
        dplyr::mutate(Subject = subject) %>% 
        dplyr::select(Subject, comparison.label, p_value) %>%
        dplyr::bind_rows(all_subjects.pvalue_table.4months_w_E)
}

#Cluster data by subjects and plot dendrogram
postscript("Clustering_dendro.LM_pvalues.WithEnergy.4months.eps",  width = 480, height = 480)
PlotDendro(all_subjects.pvalue_table.4months_w_E,
           "Cluster subjects by LM-fit p-value using 4 months of data (with energy)")
dev.off()


# Both Visits by Subject all variables 
all_subjects.pvalue_table.bothVisits = data.frame()

for(subject in c("HCR001","HCR003","HCR004","HCR006","HCR008","HCR009")) {
    subject.set <- subset(Full.with.energy, Subject == subject)
    # postscript(paste0(subject, "_LM_pvalues_BothVisits.eps"),  width = 480, height = 480)
    figure_plot = 
        PlotR2s(subject.set[, 6:dim(subject.set)[2] - 2],
                paste0("Model fit p-values for ", subject),
                "energy", color_by_pvalue = TRUE)
    # print(figure_plot)
    # dev.off()
    
    #Extract table of p-values for current subject
    subject.pvalue_table = figure_plot$data
    
    #Format and add to existing table of p-values for all subjects
    all_subjects.pvalue_table.bothVisits <-
        FormatPvalueTableForClustering(subject.pvalue_table) %>% 
        dplyr::mutate(Subject = subject) %>% 
        dplyr::select(Subject, comparison.label, p_value) %>%
        dplyr::bind_rows(all_subjects.pvalue_table.bothVisits)
}

#Cluster data by subjects and plot dendrogram
postscript("Clustering_dendro.LM_pvalues.all_measurements.Both_visits.eps",  width = 480, height = 480)
PlotDendro(all_subjects.pvalue_table.bothVisits,
           "Cluster subjects by LM-fit p-value using both visits (all measurments)")
dev.off()


#Map Subject ID to Days cutoffs for Visit 1 and Visit 2
subject_to_session.by_Days =
    data.frame(Subject = c("HCR001","HCR003","HCR004","HCR006","HCR008","HCR009"),
               Visit1 = c(43, 43, 45, 44, 45, 51),
               Visit2 = c(55, 55, 55, 57, 55, 62))

# Visit 1 by Subject all variables
all_subjects.pvalue_table.visit1 = data.frame()

for(subject in c("HCR001","HCR003","HCR004","HCR006","HCR008","HCR009")) {
    subject.set <- subset(subset(Full.with.energy, Subject == subject),
                          Days <= subset(subject_to_session.by_Days, Subject == subject)$Visit1)
    # postscript(paste0(subject, "_LM_pvalues_visit1.eps"),  width = 480, height = 480)
    figure_plot = 
        PlotR2s(subject.set[, 6:dim(subject.set)[2] - 2],
                paste0("Model fit p-values for ", subject, " (visit 1)"), "energy",
                color_by_pvalue = TRUE)
    # print(figure_plot)
    # dev.off()
    
    #Extract table of p-values for current subject
    subject.pvalue_table = figure_plot$data
    
    #Format and add to existing table of p-values for all subjects
    all_subjects.pvalue_table.visit1 <-
        FormatPvalueTableForClustering(subject.pvalue_table) %>% 
        dplyr::mutate(Subject = subject) %>% 
        dplyr::select(Subject, comparison.label, p_value) %>%
        dplyr::bind_rows(all_subjects.pvalue_table.visit1)
}

#Cluster data by subjects and plot dendrogram
postscript("Clustering_dendro.LM_pvalues.all_measurements.Visit1.eps",  width = 480, height = 480)
PlotDendro(all_subjects.pvalue_table.visit1,
           "Cluster subjects by LM-fit p-value using Visit 1 (all measurments)")
dev.off()




# Visit 2 by Subject all variables
#Note: Skip HCR008 because it's missing a lot of data from this visit
all_subjects.pvalue_table.visit2 = data.frame()

for(subject in c("HCR001","HCR003","HCR004","HCR006","HCR009")) {
    subject.set <- subset(subset(Full.with.energy, Subject == subject),
                          Days >= subset(subject_to_session.by_Days, Subject == subject)$Visit2)
    # postscript(paste0(subject, "_LM_pvalues_visit2.eps"),  width = 480, height = 480)
    figure_plot = 
        PlotR2s(subject.set[, 6:dim(subject.set)[2] - 2],
                paste0("Model fit p-values for ", subject, " (visit 2)"), "energy",
                color_by_pvalue = TRUE)
    # print(figure_plot)
    # dev.off()
    
    #Extract table of p-values for current subject
    subject.pvalue_table = figure_plot$data
    
    #Format and add to existing table of p-values for all subjects
    all_subjects.pvalue_table.visit2 <-
        FormatPvalueTableForClustering(subject.pvalue_table) %>% 
        dplyr::mutate(Subject = subject) %>% 
        dplyr::select(Subject, comparison.label, p_value) %>%
        dplyr::bind_rows(all_subjects.pvalue_table.visit2)
}

#Cluster data by subjects and plot dendrogram
postscript("Clustering_dendro.LM_pvalues.all_measurements.Visit2.eps",  width = 480, height = 480)
PlotDendro(all_subjects.pvalue_table.visit2,
           "Cluster subjects by LM-fit p-value using Visit 2 (all measurments)")
dev.off()


#Compare LM-fit p-values for each subject over both visits, visit1, visit2,
#4 months of data with energy, and 4 months of data without energy. This will
#Indicate whether the individual subjects, or the measurment periods/data types
#are having a greater affect on the cluster.

#Merge data from different measurement periods for clustering
merged_data_for_clustering <-
    dplyr::bind_rows(
        dplyr::mutate(all_subjects.pvalue_table.4months,
                      Subject = paste0(Subject, ".4months")),
        dplyr::mutate(all_subjects.pvalue_table.4months_w_E,
                      Subject = paste0(Subject, ".4months_w_E")),
        dplyr::mutate(all_subjects.pvalue_table.bothVisits,
                      Subject = paste0(Subject, ".bothVisits")),
        dplyr::mutate(all_subjects.pvalue_table.visit1,
                      Subject = paste0(Subject, ".visit1")),
        dplyr::mutate(all_subjects.pvalue_table.visit2,
                      Subject = paste0(Subject, ".visit2")))

#Cluster and plot merged data
postscript("Clustering_dendro.LM_pvalues.all_measurements.all_periods.eps",  width = 480, height = 480)    
PlotDendro(merged_data_for_clustering,
           "Cluster subjects by LM-fit p-value using all measurement periods")
dev.off()

PlotDendro.w_ggdendro(merged_data_for_clustering,
                      "Cluster subjects by LM-fit p-value using all measurement periods")





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

