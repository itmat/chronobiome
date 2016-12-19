##########################
# Amy Campbell 2016
# Plotting variability explained between each
##########################

library("ggplot2")
library("Hmisc")
library("reshape2")

###########
# Load Data
###########

aggregate_communication <-
  read.csv(
    "chronobiome/Carsten_GingerIO/
    Communication_data/Comm_data_w_results_over_4_months
    .csv"
  )
activity <-
  read.csv(
    "chronobiome/Carsten_GingerIO/Activity_data/
    Activity_data_w_results_over_4_months_1min-res.csv"
  )

# Remove missing data points(after removing Unreturned.Calls column, which in this file has huge chunks of missing values)
aggregate_communication$Unreturned.Calls <- NULL
communication_noNA <- na.omit(aggregate_communication)

View(Activity)

# Defining Functions
convert_time <- function(x) {
  #Converts time to 0-24 hours from 00:00
  if (grepl("PM", toString(x))) {
    return (as.numeric(substr(x, 11, 13)) + 12)
  }
  else{
    return(as.numeric(substr(x, 11, 13)))
  }
}

convert_date <- function(y, ord) {
  # Converts dates in dataframe y to
  # days since day 0 (10/21/2014)
  if (ord == "month_date_year") {
    return(as.Date(y, "%m/%d/%Y") - as.Date("2014-10-21"))
  }
  else{
    return(as.Date(y) - as.Date("2014-10-21"))
  }
}

Get_R2 <- function(V1, V2) {
  # Generates a simple linear regression model of V1~V2, and returns
  # R-squared for that model to represent proportion of V1 explained by V2
  V1name = colnames(V1)[1]
  V2name = colnames(V2)[1]
  mod <- lm(as.numeric(V1) ~ as.numeric(V2))
  print(V1name)
  print(V2name)
  print(mod$r.squared)
  return(c(V1name, V2name, (summary(mod)$r.squared)))
  #return(0)
}

Format_PCA_Dates <- function(v, h) {
  # Takes in visit v, hour h, and translates 'Visit' to Time Index ("hours from 10-21-2014")
  print(v)
  if (v == 1) {
    return((as.Date("2014-12-01") - as.Date("2014-10-21")) * 24 + h)
  }
  else{
    return((as.Date("2014-12-15") - as.Date("2014-10-21")) * 24 + h)
  }
}

get_Befores <- function(index) {
  #print(index)
  timeindex <- strsplit(toString(index), split = "_")
  substring <- unlist(timeindex)[2]
  timeindex <- (unlist(timeindex)[1])
  timelist <- c(seq(as.numeric(timeindex) - 5, as.numeric(timeindex)))
  timelist <- paste(timelist, substring, sep = "_")
  #print(typeof(timelist[1]))
  return(timelist)
}

Average_Befores <- function(dataframe) {
  dataframe <- dataframe
  AverageDF <- c()

  for (t in unique(dataframe$TimeSubjectIndex)) {
    beforeslist <- get_Befores(t)
    df <-
      subset(Communication_Activity, TimeSubjectIndex %in% beforeslist)
    if (dim(df)[1] == 0) {
      rowlist <- c(t, rep("NA", times = 16))
    }
    else{
      rowlist <- c(
        t,
        mean(as.numeric(df$log.Steps)),
        mean(as.numeric(df$log.Axis1)),
        mean(as.numeric(df$log.Axis2)),
        mean(as.numeric(df$log.Axis3)),
        mean(as.numeric(df$log.Lux)),
        mean(as.numeric(df$Vector.Magnitude)),
        mean(as.numeric(df$log.signal.act)),
        mean(as.numeric(df$circadian.signal.act)),
        mean(as.numeric(df$log.Mobility)),
        mean(as.numeric(df$log.Mobility.Radius)),
        mean(as.numeric(df$sqrt.Interaction.Diversity)),
        mean(as.numeric(df$SMS.Count)),
        mean(as.numeric(df$SMS.Length)),
        mean(as.numeric(df$Call.Count)),
        mean(as.numeric(df$log.signal.com)),
        mean(as.numeric(df$circadian.signal.com))
      )
      AverageDF <- rbind(AverageDF, rowlist)
    }

  }
  colnames(AverageDF) = c(
    "TimeSubjectIndex",
    "log.Steps.before",
    "log.Axis1.before",
    "log.Axis2.before",
    "log.Axis3.before",
    "log.Lux.before",
    "Vector.Magnitude.before",
    "log.signal.act.before",
    "circadian.signal.act.before",
    "log.Mobility.before",
    "log.Mobility.Radius.before",
    "sqrt.Interaction.Diversity.before",
    "SMS.Count.before",
    "SMS.Length.before",
    "Call.Count.before",
    "log.signal.com.before",
    "circadian.signal.com.before"
  )
  #View(AverageDF)
  dataframe <- merge(dataframe, AverageDF, by = "TimeSubjectIndex")
  return(dataframe)
}

get_Afters <- function(index) {
  #print(index)
  timeindex <- strsplit(toString(index), split = "_")
  substring <- unlist(timeindex)[2]
  timeindex <- (unlist(timeindex)[1])
  timelist <-
    c(seq(as.numeric(timeindex), as.numeric(timeindex) + 5))
  timelist <- paste(timelist, substring, sep = "_")
  #print(typeof(timelist[1]))
  return(timelist)
}

Average_Afters <- function(dataframe) {
  dataframe <- dataframe
  AverageDF <- c()
  for (t in unique(dataframe$TimeSubjectIndex)) {
    afterslist <- get_Afters(t)
    df <-
      subset(Communication_Activity, TimeSubjectIndex %in% afterslist)
    if (dim(df)[1] == 0) {
      rowlist <- c(t, rep("NA", times = 16))
    }
    else{
      rowlist <- c(
        t,
        mean(as.numeric(df$log.Steps)),
        mean(as.numeric(df$log.Axis1)),
        mean(as.numeric(df$log.Axis2)),
        mean(as.numeric(df$log.Axis3)),
        mean(as.numeric(df$log.Lux)),
        mean(as.numeric(df$Vector.Magnitude)),
        mean(as.numeric(df$log.signal.act)),
        mean(as.numeric(df$circadian.signal.act)),
        mean(as.numeric(df$log.Mobility)),
        mean(as.numeric(df$log.Mobility.Radius)),
        mean(as.numeric(df$sqrt.Interaction.Diversity)),
        mean(as.numeric(df$SMS.Count)),
        mean(as.numeric(df$SMS.Length)),
        mean(as.numeric(df$Call.Count)),
        mean(as.numeric(df$log.signal.com)),
        mean(as.numeric(df$circadian.signal.com))
      )
      AverageDF <- rbind(AverageDF, rowlist)
    }

  }
  colnames(AverageDF) = c(
    "TimeSubjectIndex",
    "log.Steps.after",
    "log.Axis1.after",
    "log.Axis2.after",
    "log.Axis3.after",
    "log.Lux.after",
    "Vector.Magnitude.after",
    "log.signal.act.after",
    "circadian.signal.act.after",
    "log.Mobility.after",
    "log.Mobility.Radius.after",
    "sqrt.Interaction.Diversity.after",
    "SMS.Count.after",
    "SMS.Length.after",
    "Call.Count.after",
    "log.signal.com.after",
    "circadian.signal.com.after"
  )
  #View(AverageDF)
  dataframe <- merge(dataframe, AverageDF, by = "TimeSubjectIndex")
  return(dataframe)
}

get_both <- function(index) {
  #print(index)
  timeindex <- strsplit(toString(index), split = "_")
  substring <- unlist(timeindex)[2]
  timeindex <- (unlist(timeindex)[1])
  timelist <-
    c(seq(as.numeric(timeindex) - 5 , as.numeric(timeindex) + 5))
  timelist <- paste(timelist, substring, sep = "_")
  #print(typeof(timelist[1]))
  return(timelist)
}

Average_both <- function(dataframe) {
  dataframe <- dataframe
  AverageDF <- c()

  for (t in unique(dataframe$TimeSubjectIndex)) {
    bothlist <- get_both(t)
    df <-
      subset(Communication_Activity, TimeSubjectIndex %in% bothlist)
    if (dim(df)[1] == 0) {
      rowlist <- c(t, rep("NA", times = 16))
    }
    else{
      rowlist <- c(
        t,
        mean(as.numeric(df$log.Steps)),
        mean(as.numeric(df$log.Axis1)),
        mean(as.numeric(df$log.Axis2)),
        mean(as.numeric(df$log.Axis3)),
        mean(as.numeric(df$log.Lux)),
        mean(as.numeric(df$Vector.Magnitude)),
        mean(as.numeric(df$log.signal.act)),
        mean(as.numeric(df$circadian.signal.act)),
        mean(as.numeric(df$log.Mobility)),
        mean(as.numeric(df$log.Mobility.Radius)),
        mean(as.numeric(df$sqrt.Interaction.Diversity)),
        mean(as.numeric(df$SMS.Count)),
        mean(as.numeric(df$SMS.Length)),
        mean(as.numeric(df$Call.Count)),
        mean(as.numeric(df$log.signal.com)),
        mean(as.numeric(df$circadian.signal.com))
      )
      AverageDF <- rbind(AverageDF, rowlist)
    }

  }
  colnames(AverageDF) = c(
    "TimeSubjectIndex",
    "log.Steps.both",
    "log.Axis1.both",
    "log.Axis2.both",
    "log.Axis3.both",
    "log.Lux.both",
    "Vector.Magnitude.both",
    "log.signal.act.both",
    "circadian.signal.act.both",
    "log.Mobility.both",
    "log.Mobility.Radius.both",
    "sqrt.Interaction.Diversity.both",
    "SMS.Count.both",
    "SMS.Length.both",
    "Call.Count.both",
    "log.signal.com.both",
    "circadian.signal.com.both"
  )
  #View(AverageDF)
  dataframe <- merge(dataframe, AverageDF, by = "TimeSubjectIndex")
  return(dataframe)
}



Plot <- function(dataframe) {
  r2Matrix <- c()
  for (a in colnames(dataframe)) {
    for (b in colnames(dataframe)) {
      r2Matrix <-
        rbind(r2Matrix, (Get_R2(
          as.matrix(dataframe[a]), as.matrix(dataframe[b])
        )))
    }
  }
  Plot_Matrix <- data.frame(r2Matrix)
  Plot_Matrix$X3 <- as.numeric(as.character(Plot_Matrix$X3))
  plot <-
    ggplot(data = Plot_Matrix, aes(x = X1, y = X2, fill = X3)) + geom_tile(color =
                                                                             "black") + scale_fill_gradient(low = "white", high = "steelblue") + theme(axis.text.x = element_text(angle = 90, size =
                                                                                                                                                                                    12),
                                                                                                                                                       axis.text.y = element_text(size = 13))
  return(plot)
}


# Format date and time for communication dataset
Communication_NoNA$Date <- substr(Communication_NoNA$start, 0, 10)
Communication_NoNA$Time <- convert_time(Communication_NoNA$start)

Communication_NoNA$Days <- convert_date(Communication_NoNA$Date, "")
Activity$Days <- convert_date(Activity$Date, "month_date_year")

# Set common TimeIndex and TimeSubjectIndex to synchronize the two datasets
Communication_NoNA$TimeIndex <-
  Communication_NoNA$Days * 24 + Communication_NoNA$Time
Communication_NoNA$TimeSubjectIndex <-
  paste(Communication_NoNA$TimeIndex,
        Communication_NoNA$user_id,
        sep = "_")

Activity$TimeIndex <- Activity$Days * 24 + Activity$Hour
Activity$TimeSubjectIndex <-
  paste(Activity$TimeIndex, Activity$user_id, sep = "_")

# Subset both datasets to include only hours with observations in both
Communication_Subset = subset(
  Communication_NoNA,
  Communication_NoNA$TimeSubjectIndex %in% intersect(
    Activity$TimeSubjectIndex,
    Communication_NoNA$TimeSubjectIndex
  )
)

Activity_Subset = subset(
  Activity,
  Activity$TimeSubjectIndex %in% intersect(
    Activity$TimeSubjectIndex,
    Communication_NoNA$TimeSubjectIndex
  )
)

#Join activity and communication dataframes by averaging by the hour, matching observations by TimeSubjectIndex
Communication_Activity = c()
for (x in unique(Activity_Subset$TimeSubjectIndex)) {
  sub <-
    subset(Activity_Subset, Activity_Subset$TimeSubjectIndex == x)
  subcom <-
    subset(Communication_Subset,
           Communication_Subset$TimeSubjectIndex == x)

  # Construct row to add to Communication_Activity dataframe with necessary transformations for normality)
  rowlist <- c(
    sub$TimeSubjectIndex[1],
    mean(sub$Steps),
    mean(sub$Axis1),
    mean(sub$Axis2),
    mean(sub$Axis3),
    mean(sub$Lux),
    mean(sub$Vector.Magnitude),
    mean(sub$signal),
    mean(sub$circadian.signal),
    mean(subcom$Mobility),
    mean(subcom$Mobility.Radius),
    mean(subcom$Interaction.Diversity),
    mean(subcom$SMS.Count),
    mean(subcom$SMS.Length),
    mean(subcom$Call.Count),
    mean(subcom$signal),
    mean(subcom$circadian.signal)
  )
  Communication_Activity <- rbind(Communication_Activity, rowlist)
}


colnames(Communication_Activity) <-
  c(
    "TimeSubjectIndex",
    "Steps",
    "Axis1",
    "Axis2",
    "Axis3",
    "Lux",
    "Magnitude",
    "signal.act",
    "circadian.signal.act",
    "Mobility",
    "Mobility.Radius",
    "Interaction.Diversity",
    "SMS.Count",
    "SMS.Length",
    "Call.Count",
    "signal.com",
    "circadian.signal.com"
  )
# Name columns in Transformed_Communication_Activity (log.X refers to log(X+1))

rownames(Communication_Activity) <- 1:nrow(Communication_Activity)
Communication_Activity <- data.frame(Communication_Activity)
Transformed_Communication_Activity <-
  data.frame(id = integer(dim(Communication_Activity)[1]))
View(Transformed_Communication_Activity)
################################################
View(Transformed_Communication_Activity)
Transformed_Communication_Activity["TimeSubjectIndex"] <-
  Communication_Activity$TimeSubjectIndex
Transformed_Communication_Activity["log.Steps"] <-
  log(as.numeric(Communication_Activity$Steps) + 1)
Transformed_Communication_Activity["log.Axis1"] <-
  log(as.numeric(Communication_Activity$Axis1) + 1)
Transformed_Communication_Activity["log.Axis2"] <-
  log(as.numeric(Communication_Activity$Axis2) + 1)
Transformed_Communication_Activity["log.Axis3"] <-
  log(as.numeric(Communication_Activity$Axis3) + 1)
Transformed_Communication_Activity["log.Lux"] <-
  log(as.numeric(Communication_Activity$Lux) + 1)
Transformed_Communication_Activity["Vector.Magnitude"] <-
  (Communication_Activity$Vector.Magnitude)
Transformed_Communication_Activity["log.signal.act"] <-
  log(as.numeric(Communication_Activity$signal.act) + 1)
Transformed_Communication_Activity["circadian.signal.act"] <-
  as.numeric(Communication_Activity$circadian.signal.act)
Transformed_Communication_Activity["log.Mobility"] <-
  log(as.numeric(Communication_Activity$Mobility) + 1)
Transformed_Communication_Activity["log.Mobility.Radius"] <-
  log(as.numeric(Communication_Activity$Mobility.Radius) + 1)
Transformed_Communication_Activity["sqrt.Interaction.Diversity"] <-
  sqrt(as.numeric(Communication_Activity$Interaction.Diversity))
Transformed_Communication_Activity["SMS.Count"] <-
  Communication_Activity$SMS.Count
Transformed_Communication_Activity["SMS.Length"] <-
  Communication_Activity$SMS.Length
Transformed_Communication_Activity["Call.Count"] <-
  Communication_Activity$Call.Count
Transformed_Communication_Activity["log.signal.com"] <-
  Communication_Activity$signal.com
Transformed_Communication_Activity["circadian.signal.com"] <-
  Communication_Activity$circadian.signal.com
Transformed_Communication_Activity["id"] = NULL
View(Transformed_Communication_Activity)
r2_Large_Matrix <- c()
Transformed_Communication_Activity$TimeSubjectIndex <- NULL


Plot(Transformed_Communication_Activity)
