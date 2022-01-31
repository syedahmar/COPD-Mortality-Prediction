##### Metadata ####

## File Description ##
#This file converts the data into counting process format
## Paper Title ##
# Development and Validation of a Multivariable Mortality Risk Prediction Model for COPD in Primary Care
# Reference: [TBA]
###################################### 
#Created by Syed Ahmar Shah
# Last Update: January 31, 2022
# Shared under Creative Commons Licence 
#######################################

library(survival)
# switch event values
I0=which(data$eventStatus==0)
I1=which(data$eventStatus==1)
data$eventStatus[I0]=1
data$eventStatus[I1]=0

cut_times=data$eventTime[data$eventStatus==1]

# let us divide the data into chunks before surv splitting
split_size=10000
num_parts=ceiling(length(data$patid)/split_size)

for (kpart in (1:num_parts)){
  print(paste('Processing file number:', kpart))
  index_start=1+(kpart-1)*split_size
  index_end=min( ((kpart)*split_size),length(data$patid))
  data_sel=data[index_start:index_end,]
  data_cp=survSplit(data_sel,cut=cut_times,end="eventTime",event="eventStatus",start="start")
  #data_cp_name=paste('data_cp_',kpart,sep='')
  #assign(data_cp_name,data_cp)
  save(data_cp,file=paste("data_cp_",kpart,'Rdata',sep=''))
  rm(data_sel)
  rm(data_cp)
}