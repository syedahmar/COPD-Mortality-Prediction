##### Metadata ####

## File Description ##
#Now that the data is converted into counting process format, let us process the data so that the covariates are correctly coded 
## Paper Title ##
# Development and Validation of a Multivariable Mortality Risk Prediction Model for COPD in Primary Care
# Reference: [TBA]
###################################### 
#Created by Syed Ahmar Shah
# Last Update: January 31, 2022
# Shared under Creative Commons Licence 
#######################################

library(dplyr)
library(lubridate)
chunk_size=500000;
fileno=0;
for (kfile in c(1:4)){
  load(file=paste("data_cp_",kfile,'Rdata',sep=''))
  num_chunks=ceiling(length(data_cp$patid)/chunk_size)
  for (kp in c(1:num_chunks))
  {
    print(paste('processing chunk number:',kp))
    i_start=(kp-1)*chunk_size+1;
    i_end=max( ((kp)*chunk_size),length(data_cp$patid) );
    data_cp_subset=data_cp[i_start:i_end,]
    data_cp_subset <- data_cp_subset%>%
      select(patid,indexdate,pracid,gender,birthyear,CSdeath,smokingdate,asthmadate,ARIdate,Respsympdate,Occupdate,Famhxdate,Pulmtubdate,Exerncconedate,Exercdate,IHDdate,Hyptensdate,Heartfdate,Anxietydate,Depressiondate,ses,allcausedeath,start,eventTime,eventStatus) %>%
      mutate(age_interval=round( (as.numeric(dmy(indexdate)-dmy(paste('1-7-',birthyear,sep='')),unit="days")+eventTime)/365)) %>%
      mutate(hypertension_present=ifelse( (as.numeric(dmy(Hyptensdate)-dmy(indexdate),units="days"))>eventTime | is.na(as.numeric(dmy(Hyptensdate)-dmy(indexdate),units="days")),0,1 ))%>%
      mutate(asthma_present=ifelse( (as.numeric(dmy(asthmadate)-dmy(indexdate),units="days"))>eventTime | is.na(as.numeric(dmy(asthmadate)-dmy(indexdate),units="days")),0,1 ))%>%
      mutate(depression_present=ifelse( (as.numeric(dmy(Depressiondate)-dmy(indexdate),units="days"))>eventTime | is.na(as.numeric(dmy(Depressiondate)-dmy(indexdate),units="days")),0,1 ))%>%
      mutate(HeartFailure_present=ifelse( (as.numeric(dmy(Heartfdate)-dmy(indexdate),units="days"))>eventTime | is.na(as.numeric(dmy(Heartfdate)-dmy(indexdate),units="days")),0,1 ))%>%
      mutate(IHD_present=ifelse( (as.numeric(dmy(IHDdate)-dmy(indexdate),units="days"))>eventTime | is.na(as.numeric(dmy(IHDdate)-dmy(indexdate),units="days")),0,1 ))%>%
      mutate(Anxiety_present=ifelse( (as.numeric(dmy(Anxietydate)-dmy(indexdate),units="days"))>eventTime | is.na(as.numeric(dmy(Anxietydate)-dmy(indexdate),units="days")),0,1 ))%>%
      mutate(Smoking_present=ifelse( (as.numeric(dmy(smokingdate)-dmy(indexdate),units="days"))>eventTime | is.na(as.numeric(dmy(smokingdate)-dmy(indexdate),units="days")),0,1 ))%>%
      mutate(ARI_present=ifelse( (as.numeric(dmy(ARIdate)-dmy(indexdate),units="days"))>eventTime | is.na(as.numeric(dmy(ARIdate)-dmy(indexdate),units="days")),0,1 ))%>%
      mutate(Occup_present=ifelse( (as.numeric(dmy(Occupdate)-dmy(indexdate),units="days"))>eventTime | is.na(as.numeric(dmy(Occupdate)-dmy(indexdate),units="days")),0,1 ))%>%
      mutate(Famhx_present=ifelse( (as.numeric(dmy(Famhxdate)-dmy(indexdate),units="days"))>eventTime | is.na(as.numeric(dmy(Famhxdate)-dmy(indexdate),units="days")),0,1 ))%>%
      mutate(Pulmtub_present=ifelse( (as.numeric(dmy(Pulmtubdate)-dmy(indexdate),units="days"))>eventTime | is.na(as.numeric(dmy(Pulmtubdate)-dmy(indexdate),units="days")),0,1 ))%>%
      mutate(Exernccone_present=ifelse( (as.numeric(dmy(Exerncconedate)-dmy(indexdate),units="days"))>eventTime | is.na(as.numeric(dmy(Exerncconedate)-dmy(indexdate),units="days")),0,1 ))%>%
      mutate(Exerc_present=ifelse( (as.numeric(dmy(Exercdate)-dmy(indexdate),units="days"))>eventTime | is.na(as.numeric(dmy(Exercdate)-dmy(indexdate),units="days")),0,1 ))
    
    save(data_cp_subset,file=paste('data_chunked_',(kp+fileno),'.Rdata',sep = ''))
    rm(data_cp_subset)
    gc()
  }
  rm(data_cp)
  gc()
  fileno=fileno+num_chunks
}