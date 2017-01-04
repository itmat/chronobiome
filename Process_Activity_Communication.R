######################################################
# Amy Campbell 2016
# Plotting variability explained between each variable
######################################################

library("Hmisc")
library("reshape2")
library("dplyr")
library("readr")
library("chron")

###############
# Load Datasets
###############

aggregate.communication <-
  read.csv(
    "chronobiome/Carsten_GingerIO/Communication_data/Comm_data_w_results_over_4_months.csv"
  )
aggregate.communication$Unreturned.Calls <- NULL
communication.NoNA <- na.omit(aggregate.communication)

activity <-
  read.csv("chronobiome/Carsten_GingerIO/Activity_data/Activity_data_w_results_over_4_months_1min-res.csv")

#######################################################
# Process data, combine activity and communication sets
#######################################################

# Standardize dates and times to be "hours from 0th hour of 2014-10-21"
communication.NoNA$Date <- substr(communication.NoNA$start, 0, 10)
communication.NoNA$Days <- chron(communication.NoNA$Date, 
                                 format=c(dates="y-m-d")) - chron("2014-10-21", 
                                 format=c(dates="y-m-d"))
activity$Days <- chron(as.character(activity$Date), 
                       format=c(dates="m/d/y")) - chron("2014-10-21", 
                       format=c(dates="y-m-d"))

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
# also present in the activity dataset, take average of each variable at 
# each Subject/Hour combination
communication.df <- communication.NoNA %>% 
  dplyr::filter(TimeSubjectIndex %in%
                  intersect(activity$TimeSubjectIndex,
                            communication.NoNA$TimeSubjectIndex)) %>%
  dplyr::group_by(TimeSubjectIndex) %>%
  dplyr::summarise(Mobility = mean(Mobility),
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

# Subset activity dataset to include only TimeSubjectIndex recordings
# also present in the communication dataset, take average of each variable at 
# each Subject/Hour combination
activity.df <- activity %>% 
 dplyr::filter(TimeSubjectIndex %in%
               intersect(activity$TimeSubjectIndex,communication.NoNA$TimeSubjectIndex)) %>%
  dplyr::group_by(TimeSubjectIndex) %>%
  dplyr::summarise(Steps = mean(Steps), Axis1 = mean(Axis1),
                   Axis2 = mean(Axis2), Axis3 = mean(Axis3),
                   Luminosity = mean(Lux), activity.Vector.Magnitude = mean(Vector.Magnitude),
                   signal.act = mean(signal),
                   circadian.signal.act = mean(circadian.signal),
                   activity.amplitude=mean(instantaneous.amplitude),
                   activity.period = mean(instantaneous.period),
                   activity.phase = mean(instantaneous.phase))

communication.activity <- dplyr::full_join(activity.df, communication.df, by="TimeSubjectIndex")

rownames(communication.activity) <- 1:nrow(communication.activity)
communication.activity <- data.frame(communication.activity)


# Transform variables in communication.activity to improve normality
# In transformed.data.activity, log.X refers to transformation log(X+1)
transformed.communication.activity <- communication.activity %>%
  dplyr::group_by(TimeSubjectIndex) %>%
  dplyr::summarise(log.Steps = log(Steps + 1),
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

# Save dataframes as .csv files

write.table(communication.activity, "communication.activity.csv", sep=",")
write.table(transformed.communication.activity, "transformed.communication.activity.csv", sep = ",")


