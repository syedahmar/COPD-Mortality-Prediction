##### Metadata ####

## File Description ##
#This file loads the data and pre-processes it to get the death timings right
## Paper Title ##
# Development and Validation of a Multivariable Mortality Risk Prediction Model for COPD in Primary Care
# Reference: [TBA]
###################################### 
#Created by Syed Ahmar Shah
# Last Update: January 31, 2022
# Shared under Creative Commons Licence 
#######################################

library(lubridate)
data = read.csv(file="datasets/UpdatedData.csv")
# let us now create a survival time column
deathDays=dmy(data$deathdate)-dmy(data$indexdate);
deathDays=as.numeric(deathDays,units="days")
data$eventTime=deathDays
time_10years=10*365
Index_notDied=which(is.na(deathDays))
Index_DiedAfter10years= which(deathDays>time_10years)
data$eventTime[Index_DiedAfter10years]=time_10years # we follow up until 10 years
data$eventTime[Index_notDied]=dmy(data$regend[Index_notDied])-dmy(data$indexdate[Index_notDied])
data$eventStatus = matrix(0,nrow(data),1)
data$eventStatus[c(Index_notDied,Index_DiedAfter10years)]=1
Index_FollowUpGreaterthan10yrs=which(data$eventTime>time_10years)
data$eventTime[Index_FollowUpGreaterthan10yrs]=time_10years