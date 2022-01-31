##### Metadata ####

## File Description ##
# Fit the Cox Regression Model and internally validate it with a processing pipeline
## Paper Title ##
# Development and Validation of a Multivariable Mortality Risk Prediction Model for COPD in Primary Care
# Reference: [TBA]
###################################### 
#Created by Syed Ahmar Shah
# Last Update: January 31, 2022
# Shared under Creative Commons Licence 
#######################################

rm(list=ls())
library(dplyr)
library(rBayesianOptimization)
library(survival)
library(SuppDists)
load('data_all_v5.Rdata')

# convert the new additional data into factors
data_all$exacerbations=as.factor(data_all$exacerbations)
data_all$mMRCScale=as.factor(data_all$mMRCScale)
data_all$breathless=as.factor(data_all$breathless)
data_all$hospitalised=as.factor(data_all$hospitalised)

data_all$smoking2=relevel(data_all$smoking2,ref="never")
data_all$exacerbations=relevel(data_all$exacerbations,ref="0")
data_all$mMRCScale=relevel(data_all$mMRCScale,ref="0")
data_all$breathless=relevel(data_all$breathless,ref="0")
data_all$hospitalised=relevel(data_all$hospitalised,ref="0")

PatientIDs = unique(data_all$patid)
num_subcohorts=8
num_kfolds=10
set.seed(23) # for reproducibility
patid_subcohorts_indices=KFold(PatientIDs, nfolds = num_subcohorts, stratified = FALSE, seed = 0)
#this loop will iterate through each subsample of cohort (each subsample is a random sampling from the original cohort)
# initialize an empty dataframe where the key results will be stored
validation_overall <- data.frame(cohort=character(),
                                 fold=character(),
                                 patients_train=character(),
                                 patients_test=character(),
                                 bestfit_train=character(),
                                 bestfit_test=character(),
                                 cindex_train=character(),
                                 cindex_test=character(),
                                 dindex_train=character(),
                                 dindex_test=character())
# create a new directory where results will be stored
ResultsDirectory='Results_Validation'
dir.create(ResultsDirectory)
for (i_cohort in c(1:num_subcohorts)){
  print(paste('Processing cohort number: ',i_cohort))
  
  # get the patient ids of the subsample now
  PatientIDs_Cohort=PatientIDs[unlist(patid_subcohorts_indices[i_cohort])]
  data_cohort <- data_all %>%
    filter(patid %in% PatientIDs_Cohort)
  # now for the selected cohort, let us do a 10-fold cross validation
  #set.seed(50) # for reproducibility
  patid_kfold_indices=KFold(PatientIDs_Cohort, nfolds = num_kfolds, stratified = FALSE, seed = 0)
  for (i_fold in c(1:num_kfolds)){
    print(paste('Processing Cohort:',i_cohort,' and Fold:',i_fold))
    # get the training and testing set now
    patient_train_IDs=PatientIDs_Cohort[unlist(patid_kfold_indices[-i_fold])]
    patient_test_IDs=PatientIDs_Cohort[unlist(patid_kfold_indices[i_fold])]
    data_train <- data_all %>%
      filter(patid %in% patient_train_IDs) 
    data_test <- data_all %>%
      filter(patid %in% patient_test_IDs)
    #remove any NA values (found in ses for some cases)
    data_test=data_test[complete.cases(data_test),]
    data_train=data_train[complete.cases(data_train),]
    # now that we have the training and testing, let us fit the cox model to the training data
    Extended_Surv_object_train = Surv(data_train$start, data_train$eventTime,data_train$eventStatus)
    coxph_extended_train = coxph(Extended_Surv_object_train~gender+age_interval+smoking2+hypertension_present+asthma_present+depression_present+HeartFailure_present+IHD_present+Anxiety_present+ARI_present+Occup_present+Famhx_present+Pulmtub_present+Exernccone_present+Exerc_present+ses+breathless+exacerbations+mMRCScale+hospitalised+cluster(patid),data=data_train)
    #now that we have run the model, it is time to validate it on the test set
    #1: Get the PI score on the test set
    PI_score_test = (coxph_extended_train$coefficients["gendermale"]*(data_test$gender=='male'))+(coxph_extended_train$coefficients['age_interval']*(data_test$age_interval))+(coxph_extended_train$coefficients['smoking2current']*(data_test$smoking2=='current'))+(coxph_extended_train$coefficients['smoking2former']*(data_test$smoking2=='former'))+(coxph_extended_train$coefficients['hypertension_present']*data_test$hypertension_present)+(coxph_extended_train$coefficients['asthma_present']*data_test$asthma_present)+(coxph_extended_train$coefficients['depression_present']*data_test$depression_present)+(coxph_extended_train$coefficients['HeartFailure_present']*data_test$HeartFailure_present)+(coxph_extended_train$coefficients['IHD_present']*data_test$IHD_present)+(coxph_extended_train$coefficients['Anxiety_present']*data_test$Anxiety_present)+(coxph_extended_train$coefficients['ARI_present']*data_test$ARI_present)+(coxph_extended_train$coefficients['Occup_present']*data_test$Occup_present)+(coxph_extended_train$coefficients['Famhx_present']*data_test$Famhx_present)+(coxph_extended_train$coefficients['Pulmtub_present']*data_test$Pulmtub_present)+(coxph_extended_train$coefficients['Exernccone_present']*data_test$Exernccone_present)+(coxph_extended_train$coefficients['Exerc_present']*data_test$Exerc_present)+(coxph_extended_train$coefficients['ses2nd_IMD_Quintile']*(data_test$ses=='2nd_IMD_Quintile'))+(coxph_extended_train$coefficients['ses3rd_IMD_Quintile']*(data_test$ses=='3rd_IMD_Quintile'))+(coxph_extended_train$coefficients['ses4th_IMD_Quintile']*(data_test$ses=='4th_IMD_Quintile'))+(coxph_extended_train$coefficients['ses5th_IMD_Quintile']*(data_test$ses=='5th_IMD_Quintile'))+(coxph_extended_train$coefficients['breathless1']*(data_test$breathless=='1'))+(coxph_extended_train$coefficients['exacerbations1']*(data_test$exacerbations=='1'))+(coxph_extended_train$coefficients['exacerbations2']*(data_test$exacerbations=='2'))+(coxph_extended_train$coefficients['mMRCScale2']*(data_test$mMRCScale=='2'))+(coxph_extended_train$coefficients['mMRCScale3']*(data_test$mMRCScale=='3'))+(coxph_extended_train$coefficients['mMRCScale8']*(data_test$mMRCScale=='8'))+(coxph_extended_train$coefficients['hospitalised1']*(data_test$hospitalised=='1'))
    
    PI_score_train = (coxph_extended_train$coefficients["gendermale"]*(data_train$gender=='male'))+(coxph_extended_train$coefficients['age_interval']*(data_train$age_interval))+(coxph_extended_train$coefficients['smoking2current']*(data_train$smoking2=='current'))+(coxph_extended_train$coefficients['smoking2former']*(data_train$smoking2=='former'))+(coxph_extended_train$coefficients['hypertension_present']*data_train$hypertension_present)+(coxph_extended_train$coefficients['asthma_present']*data_train$asthma_present)+(coxph_extended_train$coefficients['depression_present']*data_train$depression_present)+(coxph_extended_train$coefficients['HeartFailure_present']*data_train$HeartFailure_present)+(coxph_extended_train$coefficients['IHD_present']*data_train$IHD_present)+(coxph_extended_train$coefficients['Anxiety_present']*data_train$Anxiety_present)+(coxph_extended_train$coefficients['ARI_present']*data_train$ARI_present)+(coxph_extended_train$coefficients['Occup_present']*data_train$Occup_present)+(coxph_extended_train$coefficients['Famhx_present']*data_train$Famhx_present)+(coxph_extended_train$coefficients['Pulmtub_present']*data_train$Pulmtub_present)+(coxph_extended_train$coefficients['Exernccone_present']*data_train$Exernccone_present)+(coxph_extended_train$coefficients['Exerc_present']*data_train$Exerc_present)+(coxph_extended_train$coefficients['ses2nd_IMD_Quintile']*(data_train$ses=='2nd_IMD_Quintile'))+(coxph_extended_train$coefficients['ses3rd_IMD_Quintile']*(data_train$ses=='3rd_IMD_Quintile'))+(coxph_extended_train$coefficients['ses4th_IMD_Quintile']*(data_train$ses=='4th_IMD_Quintile'))+(coxph_extended_train$coefficients['ses5th_IMD_Quintile']*(data_train$ses=='5th_IMD_Quintile'))+(coxph_extended_train$coefficients['breathless1']*(data_train$breathless=='1'))+(coxph_extended_train$coefficients['exacerbations1']*(data_train$exacerbations=='1'))+(coxph_extended_train$coefficients['exacerbations2']*(data_train$exacerbations=='2'))+(coxph_extended_train$coefficients['mMRCScale2']*(data_train$mMRCScale=='2'))+(coxph_extended_train$coefficients['mMRCScale3']*(data_train$mMRCScale=='3'))+(coxph_extended_train$coefficients['mMRCScale8']*(data_train$mMRCScale=='8'))+(coxph_extended_train$coefficients['hospitalised1']*(data_train$hospitalised=='1'))
    
    ##### Validation, Method 1 (Regression on PI) ####
    #test set
    PI_score_test_centered=PI_score_test-mean(PI_score_test)
    test_dataframe_m1=data.frame(c(data_test$patid),c(PI_score_test_centered))
    names(test_dataframe_m1)=c('patid','PI_score_centered')
    Surv_object_test = Surv(data_test$start, data_test$eventTime,data_test$eventStatus)
    coxph_bestfit_test = coxph(Surv_object_test~PI_score_test_centered+cluster(patid),data=test_dataframe_m1)
    coxph_bestfit_test$coefficients
    #train set
    PI_score_train_centered=PI_score_train-mean(PI_score_train)
    train_dataframe_m1=data.frame(c(data_train$patid),c(PI_score_train_centered))
    names(train_dataframe_m1)=c('patid','PI_score_centered')
    Surv_object_train = Surv(data_train$start, data_train$eventTime,data_train$eventStatus)
    coxph_bestfit_train = coxph(Surv_object_train~PI_score_train_centered+cluster(patid),data=train_dataframe_m1)
    coxph_bestfit_train$coefficients
    
    #Validation, Method 3 (c-index)
    test_cindex=survConcordance(Surv_object_test ~PI_score_test)
    test_cindex$concordance
    train_cindex=survConcordance(Surv_object_train ~PI_score_train)
    train_cindex$concordance
    ###### Validation, Method 3 (d-index) #####
    #test set
    rankits_test=normOrder(length(PI_score_test)) #get the rankits 
    rank_test=rank(PI_score_test) # rank the PI score
    PI_rankits_test=rankits_test[rank_test] # find the associated rankits of the PI score
    PI_rankits_test_scaled=PI_rankits_test/(sqrt(8/pi)) # scale the rankits
    dindex_test_coxph <- coxph(Surv_object_test ~ PI_rankits_test_scaled)
    d.index_test=dindex_test_coxph$coefficients # assign the coefficient found from res.cox to d.index
    d_index_test_final=((d.index_test^2)/(1.596)^2)/(1.645+((d.index_test^2)/(1.596^2)) )
    #train set
    rankits_train=normOrder(length(PI_score_train)) #get the rankits 
    rank_train=rank(PI_score_train) # rank the PI score
    PI_rankits_train=rankits_train[rank_train] # find the associated rankits of the PI score
    PI_rankits_train_scaled=PI_rankits_train/(sqrt(8/pi)) # scale the rankits
    dindex_train_coxph <- coxph(Surv_object_train ~ PI_rankits_train_scaled)
    d.index_train=dindex_train_coxph$coefficients # assign the coefficient found from res.cox to d.index
    d_index_train_final=((d.index_train^2)/(1.596)^2)/(1.645+((d.index_train^2)/(1.596^2)) )
    
    new_record=data.frame(c(i_cohort),c(i_fold),c(length(unique(data_train$patid))),c(length(unique(data_test$patid))),c(coxph_bestfit_train$coefficients),c(coxph_bestfit_test$coefficients),c(train_cindex$concordance),c(test_cindex$concordance),c(d_index_train_final),c(d_index_test_final))
    names(new_record)=c("cohort","fold","patients_train","patients_test","bestfit_train","bestfit_test","cindex_train","cindex_test","dindex_train","dindex_test")
    validation_overall = rbind(validation_overall,new_record)
    save(coxph_extended_train,file=paste(ResultsDirectory,'/coxph_extended_train_cohort',i_cohort,'_fold',i_fold,'.Rdata',sep=''))
    #End of fold loop
    rm(data_train)
    rm(data_test)
    rm(coxph_extended_train)
    rm(coxph_bestfit_test)
    rm(coxph_bestfit_train)
    rm(dindex_test_coxph)
    rm(dindex_train_coxph)
    rm(test_dataframe_m1)
    rm(train_dataframe_m1)
    rm(list=setdiff(ls(), c("i_cohort","i_fold","data_all","validation_overall","PatientIDs_Cohort","patid_kfold_indices","num_kfolds","num_subcohorts","ResultsDirectory","PatientIDs","patid_subcohorts_indices")))
  }
  
  #clear the variables not required anymore to save memory  
  rm(data_cohort)
  #end of cohort loop
}
save(validation_overall,file=paste(ResultsDirectory,'validation_initial.Rdata'))

