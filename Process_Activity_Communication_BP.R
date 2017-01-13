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

# Unreturned calls has a very high occurance of NA values
aggregate.communication$Unreturned.Calls <- NULL

# Having removed Unreturned Calls, eliminate rows with missing values
communication.NoNA <- na.omit(aggregate.communication)


activity <-
  read.csv("chronobiome/Carsten_GingerIO/Activity_data/Activity_data_w_results_over_4_months_1min-res.csv")

# Blood pressure/ Heartrate data
bloodpressure.HR <- read.csv("chronobiome/Carsten_GingerIO/Biometric_data/BP_data/Raw_bp_data_for_figures.txt", sep="\t")
bloodpressure.HR <- na.omit(bloodpressure.HR)



#######################################################
# Process data, combine activity and communication sets
#######################################################
# Standardize dates and times to be "hours from 0th hour of 2014-10-21"
communication.NoNA$Date <- substr(communication.NoNA$start, 0, 10)
communication.NoNA$Days <- chron(as.character(communication.NoNA$Date), 
                                 format=c(dates="y-m-d")) - chron("2014-10-21", 
                                                                  format=c(dates="y-m-d"))
activity$Days <- chron(as.character(activity$Date), 
                       format=c(dates="m/d/y")) - chron("2014-10-21", 
                                                        format=c(dates="y-m-d"))
bloodpressure.HR$Date <- substr(bloodpressure.HR$Start.Time, 0, 10)
bloodpressure.HR$Days <- chron(as.character(bloodpressure.HR$Date), 
                               format=c(dates="y-m-d")) - chron("2014-10-21", 
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

bloodpressure.HR$TimeIndex <- bloodpressure.HR$Days * 24 + bloodpressure.HR$Hour
bloodpressure.HR$TimeSubjectIndex <-
  paste(bloodpressure.HR$TimeIndex, bloodpressure.HR$Subject, sep = "_")

heartrate <- subset(bloodpressure.HR, Data.type=="BP.HR")
diastolic <- subset(bloodpressure.HR, Data.type == "BP.Diastolic")
systolic <- subset(bloodpressure.HR, Data.type == "BP.Systolic")
arterialpressure <- subset(bloodpressure.HR, Data.type == "BP.MAP")
pulsepressure <- subset(bloodpressure.HR, Data.type == "BP.PP")


# Parse bloodpressure.HR dataset
heartrate["Data.type"] <- NULL
diastolic["Data.type"] <- NULL
systolic["Data.type"] <- NULL
arterialpressure["Data.type"] <- NULL
pulsepressure["Data.type"] <- NULL

colnames(heartrate) <- c("Subject", "Start.Time", "Month", "Day", "Day_of_Week", "Hour", "heart.rate","Date","Days", "Times", "TimeSubjectIndex")
colnames(diastolic) <- c("Subject", "Start.Time", "Month", "Day", "Day_of_Week", "Hour", "diastolic.bp", "Date","Days", "Times", "TimeSubjectIndex")
colnames(systolic) <- c("Subject", "Start.Time", "Month", "Day", "Day_of_Week", "Hour", "systolic.bp", "Date","Days", "Times", "TimeSubjectIndex")
colnames(arterialpressure) <- c("Subject", "Start.Time", "Month", "Day", "Day_of_Week", "Hour", "arterial.pressure", "Date","Days", "Times", "TimeSubjectIndex")
colnames(pulsepressure) <- c("Subject", "Start.Time", "Month", "Day", "Day_of_Week", "Hour", "pulse.pressure", "Date","Days", "Times", "TimeSubjectIndex")


heartrate.df <- heartrate %>%
  dplyr::group_by(TimeSubjectIndex) %>%
  dplyr::summarise_if(is.numeric, mean)

diastolic.df <- diastolic %>%
  dplyr::group_by(TimeSubjectIndex) %>%
  dplyr::summarise(diastolic.bp = mean(diastolic.bp))

systolic.df <- systolic %>%
  dplyr::group_by(TimeSubjectIndex) %>%
  dplyr::summarise(systolic.bp = mean(systolic.bp))

arterial.df <- arterialpressure %>%
  dplyr::group_by(TimeSubjectIndex) %>%
  dplyr::summarise(arterial.pressure = mean(arterial.pressure))

pulse.df <- pulsepressure %>%
  dplyr::group_by(TimeSubjectIndex) %>%
    dplyr::summarise(pulse.pressure = mean(pulse.pressure))

Heartrate.BP <- heartrate.df %>%
  dplyr::full_join(diastolic.df, by="TimeSubjectIndex") %>%
  dplyr::full_join(systolic.df, by="TimeSubjectIndex") %>%
  dplyr::full_join(arterial.df, by="TimeSubjectIndex") %>%
  dplyr::full_join(pulse.df, by="TimeSubjectIndex")
  
Heartrate.BP <- Heartrate.BP[c("TimeSubjectIndex", "Days","Times", "heart.rate","diastolic.bp","systolic.bp", "arterial.pressure", "pulse.pressure")]

write.csv(Heartrate.BP, "heartrate.bp.csv")

# Subset communication dataset to include only TimeSubjectIndex recordings
# also present in the activity dataset, take average of each variable at 
# each Subject/Hour combination

communication.df <- communication.NoNA %>% 
  dplyr::filter(TimeSubjectIndex %in%
                  intersect(activity$TimeSubjectIndex, communication.NoNA$TimeSubjectIndex)) %>%
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

communication.activity <- dplyr::full_join(activity.df, communication.df, by="TimeSubjectIndex")


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
write.csv(communication.activity, "communication.activity.csv")
write.csv(transformed.communication.activity, "transformed.communication.activity.csv")
