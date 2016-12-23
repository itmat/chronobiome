######################################################
# Amy Campbell 2016
# Plotting variability explained between each variable
######################################################

library("ggplot2")
library("Hmisc")
library("reshape2")
library("dplyr")
library("readr")
library("chron")

####################
# Define Functions
####################

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

PlotR2s <- function(dataframe) {
  # Generates geom_tile plot of the % variability explained (R2) for each variable by each other variable
  # :param: dataframe - Dataframe to plot (each variable in the dataframe will be compared to each other
  #  dataframe, so these must be numeric)
  r2Matrix <- c()
  for (a in colnames(dataframe)) {
    for (b in colnames(dataframe)) {
      r2Matrix <-
        rbind(r2Matrix, (GetR2(
          as.matrix(dataframe[a]), as.matrix(dataframe[b])
        )))
    }
  }
  PlotMatrix <- data.frame(r2Matrix)
  colnames(PlotMatrix) <- c("factor.1", "factor.2", "variability.explained")
  PlotMatrix$variability.explained <- as.numeric(as.character(PlotMatrix$variability.explained))
  plot <-
    ggplot(data = PlotMatrix, aes(x = factor.1, y = factor.2, fill = variability.explained)) +
    geom_tile(color ="black") + scale_fill_gradient(low="white", high="steelblue4") +
    theme(axis.text.x = element_text(angle = 90, size =11), axis.text.y = element_text(size = 11)) +
    ggtitle("Variability Explained in activity and Communication Data") 

  return(plot)
}

###############
# Load Datasets
###############

aggregate.communication <-
  read.csv(
    "chronobiome/Carsten_GingerIO/Communication_data/Comm_data_w_results_over_4_months.csv"
  )
aggregate.communication$Unreturned.Calls <- NULL
communication_noNA <- na.omit(aggregate.communication)

activity <-
  read.csv("chronobiome/Carsten_GingerIO/activity_data/activity_data_w_results_over_4_months_1min-res.csv")

# Format date and time
communication.NoNA$Date <- substr(communication.NoNA$start, 0, 10)
communication.NoNA$Days <- chron(communication.NoNA$Date, format=c(dates="y-m-d")) - chron("2014-10-21", format=c(dates="y-m-d"))
activity$Days <- chron(as.character(activity$Date), format=c(dates="m/d/y")) - chron("2014-10-21", format=c(dates="y-m-d"))

# Set common TimeIndex and TimeSubjectIndex to synchronize the two datasets
communication.NoNA$TimeIndex <-
  communication.NoNA$Days * 24 + communication.NoNA$Hour
communication.NoNA$TimeSubjectIndex <-
  paste(communication.NoNA$TimeIndex,
        communication.NoNA$user_id,
        sep = "_")

activity$TimeIndex <- activity$Days * 24 + activity$Hour
activity$TimeSubjectIndex <-
  paste(activity$TimeIndex, activity$user_id, sep = "_")

# Subset communication dataset to include only TimeSubjectIndex recordings 
# also present in the activity dataset
communication.subset = subset(
  communication.NoNA,
  communication.NoNA$TimeSubjectIndex %in% intersect(
    activity$TimeSubjectIndex,
    communication.NoNA$TimeSubjectIndex
  )
)
# Average 
communication.subset <- group_by(communication.subset, TimeSubjectIndex)
communication.df <- summarise(communication.subset,
                              Mobility = mean(Mobility),
                              Mobility.Radius = mean(Mobility.Radius),
                              Interaction.Diversity = mean(Interaction.Diversity),
                              SMS.Count = mean(SMS.Count),
                              SMS.Length = mean(SMS.Length),
                              Call.Count = mean(Call.Count),
                              signal.com= mean(signal),
                              circadian.signal.com = mean(circadian.signal),
                              Communication.amplitude = mean(instantaneous.amplitude),
                              Communication.period = mean(instantaneous.period),
                              Communication.phase = mean(instantaneous.phase))

activity_Subset <- subset(
  activity,
  activity$TimeSubjectIndex %in% intersect(
    activity$TimeSubjectIndex,
    communication.NoNA$TimeSubjectIndex
  )
)

activity_Subset <- group_by(activity_Subset, TimeSubjectIndex)
activity.df <- summarise(activity_Subset,
                   Steps = mean(Steps), Axis1 = mean(Axis1), 
                   Axis2 = mean(Axis2), Axis3 = mean(Axis3), 
                   Luminosity = mean(Lux), activity.Vector.Magnitude = mean(Vector.Magnitude),
                   signal.act = mean(signal),
                   circadian.signal.act = mean(circadian.signal),
                   activity.amplitude=mean(instantaneous.amplitude),
                   activity.period = mean(instantaneous.period),
                   activity.phase = mean(instantaneous.phase))

communication.activity <- full_join(activity.df, communication.df, by="TimeSubjectIndex")

# Name columns in transformed.communication.activity (log.X refers to log(X+1))

rownames(communication.activity) <- 1:nrow(communication.activity)
communication.activity <- data.frame(communication.activity)
write.table(communication.activity, "communication.activity.csv", sep=",")

communication.activity <- group_by(communication.activity,TimeSubjectIndex)
transformed.communication.activity <- summarise(communication.activity,
                                                log.Steps = log(Steps + 1),
                                                log.Axis1 = log(Axis1 + 1),
                                                log.Axis2 = log(Axis2 + 1),
                                                log.Axis3 = log(Axis3 + 1),
                                                log.Luminosity = log(Luminosity + 1),
                                                log.activity.Vector.Magnitude = log(activity.Vector.Magnitude +1),
                                                log.signal.act = log(signal.act + 1),
                                                circadian.signal.act = circadian.signal.act,
                                                activity.amplitude = activity.amplitude,
                                                activity.period = activity.period,
                                                activity.phase = activity.phase,
                                                log.Mobility = log(Mobility + 1),
                                                log.Mobility.Radius = log(Mobility.Radius + 1),
                                                sqrt.Interaction.Diversity = sqrt(Interaction.Diversity),
                                                SMS.Count = SMS.Count,
                                                SMS.Length = SMS.Length,
                                                Call.Count = Call.Count,
                                                log.signal.com = log(signal.com + 1),
                                                circadian.signal.com = circadian.signal.com,
                                                Communication.amplitude = Communication.amplitude,
                                                Communication.period = Communication.period, 
                                                Communication.phase = Communication.phase)

transformed.communication.activity$TimeSubjectIndex <- NULL
jpeg('CommunicationactivityPlot.jpg')
PlotR2s(transformed.communication.activity)
dev.off()

