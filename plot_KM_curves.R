##### Metadata ###################################################
#This file plots survival curves using the formated data for 
#deriving the COPD Mortality Risk Prediction Model
################################################################## 
#Created by Syed Ahmar Shah
# Last Update: January 31, 2022
# Shared under Creative Commons Licence 
##################################################################


library(survival)
library(survminer)
# Create a survival object (mMRC Score)
daysInYear=365
Surv_object_1 = Surv(data$eventTime/daysInYear,data$eventStatus==0)
Surv_fit_1=survfit(Surv_object_1~1)

kmfit_1=(survfit(Surv_object_1~data$exacerbations))
plot(kmfit_1,col=(1:4),xlab="survival time in years"
     ,ylab="survival probability")
legend("topright",rownames(summary(kmfit_1)$table),lty=c("solid","dashed"),col=(1:4))
summary(kmfit_1,times=c(0:10))
g1=ggsurvplot(kmfit_1,data=data,risk.table=TRUE,risk.table.height = 0.3,ggtheme = theme_bw(),xlab="time(years)",legend.labs=c('exacerbations=0','exacerbations=1','exacerbations\u2265 2') )
#g1=ggsurvplot(kmfit_1,data=data,risk.table=TRUE,risk.table.height = 0.3,ggtheme = theme_bw(),xlab="time(years)",legend.labs=c('mMRC dyspnoea score(missing)','0-1','2','3-4') )
g1

# Create a survival object (Exacerbations)
daysInYear=365
Surv_object_1 = Surv(data$eventTime/daysInYear,data$eventStatus==0)
Surv_fit_1=survfit(Surv_object_1~1)

kmfit_1=(survfit(Surv_object_1~data$hospitalised))
plot(kmfit_1,col=(1:4),xlab="survival time in years"
     ,ylab="survival probability")
legend("topright",rownames(summary(kmfit_1)$table),lty=c("solid","dashed"),col=(1:4))
summary(kmfit_1,times=c(0:10))
ggsurvplot(kmfit_1,data=data,risk.table=TRUE,risk.table.height = 0.3,ggtheme = theme_bw(),xlab="time(years)" )