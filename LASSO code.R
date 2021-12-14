rm(list = ls())

# multiply-imputed by chained equations
library(readxl)
data1 <- read_excel("./data/original_data1.xlsx")
library(lattice)
library(MASS)
library(mice)
imp =mice(data1,m=5)
fit=with(imp,lm(status~time,data=data1))
result=complete(imp,action=3)


# LASSO and a ten-fold cross validation
data <- read_excel("./data/data.xlsx")
library(Matrix)
library(glmnet)
library(survival)
str(data)
# Change the age to a numeric
fix(data) # Select the numeric after clicking on age
x <- as.matrix(data[,5:78])
y <- as.matrix(Surv(data$time,data$status))
lasso <- glmnet(x,y,alpha=1,family="cox")
print(lasso)
plot(lasso,xvar="lambda",label = "TRUE")
set.seed(2)
cv.lasso <- cv.glmnet(x,y,family="cox",alpha=1,nfolds=10)
plot(cv.lasso)
coef(cv.lasso,s="lambda.min")
lasso <- cv.lasso$lambda.min
lasso
log(lasso)


# Univariate Cox regression
library(survival)
library(plyr)
y<- Surv(time = data$time,event = data$status==1)
Uni_cox_model<-
  function(x){
    FML <- as.formula(paste0 ("y~",x))
    cox<- coxph(FML,data=data)
    cox1<-summary(cox)
    HR <- round(cox1$coefficients[,2],3)
    PValue <- round(cox1$coefficients[,5],3)
    CI5 <-round(cox1$conf.int[,3],3)
    CI95 <-round(cox1$conf.int[,4],3)
    Uni_cox_model<- data.frame(
      names <-rownames(cox1$conf.int),
      'HR' = HR,
      'CI5' = CI5,
      'CI95' = CI95,
      'P' = PValue)
    return(Uni_cox_model)
  }  

variable.names<- colnames(data)[c(4,6,11,21,27,28,34,37,38,44,50,56,57,60,71)];variable.names
Uni_cox <- lapply(variable.names, Uni_cox_model)
Uni_cox<- ldply(Uni_cox,data.frame)
Uni_cox$HR.CI95<-paste0(Uni_cox$HR," (",Uni_cox$CI5,'-',Uni_cox$CI95,")");Uni_cox



Multivariable  Cox regression
library(survival)
mul_cox<-coxph(Surv(time,status==1)~
                 age+RA+Immunosuppressive.agents+DLCO+RVD+RAA+PASP+LVEF+CRP+BNP+AST+GGT+ALB+Honeycombing,
               data=data
);summary(mul_cox)
cox<-summary(mul_cox) 
cox$coefficients    
cox$conf.int  
mul_HR<- round(cox$coefficients[,2],3) 
mul_PValue<- round(cox$coefficients[,5],3) 
mul_CI1<-round(cox$conf.int[,3],3)
mul_CI2<-round(cox$conf.int[,4],3)
mul_CI95<-paste(mul_CI1,'-',mul_CI2)
mul_cox1 <- data.frame("HR" =mul_HR,
                       "CI95" =mul_CI95,
                       "P"=mul_PValue);mul_cox1



# Nomogram
suppressMessages(library(rms))
dd <- datadist(data)
options(datadist="dd")
cox_nomo1 <-  cph(Surv(time,status)~age+RA+DLCO+RVD+RAA+Honeycombing+AST+ALB+Immunosuppressive.agents,data=data,x=T,y=T,surv=T)
surv <- Survival(cox_nomo1)
surv1 <- function(x)surv(1*36,lp=1-x)#"lp=1-x)
surv2 <- function(x)surv(1*60,lp=1-x)
nomo_2a<-nomogram(cox_nomo1, fun=list(surv1,surv2), lp=F,funlabel =c("3-year Death", "5-year Death"),fun.at =c(0.05, seq(0.1,0.9, by=0.1), 0.95))
plot(nomo_2a, col.grid=c("pink","cyan"),xfrac = 0.3,cex.var = 1,cex.axis = 1,lmgp = 0.3)


Harrell¡¯s C index
data$p_lp <- predict(cox_nomo1, data, type="lp")
c_harrell_nomogram <- (cph(Surv(time,status)~p_lp, data=data,x=TRUE,y=TRUE)$stats["Dxy"]+1)/2
c_harrell_nomogram

cox_GAP <-  cph(Surv(time,status)~GAP,data=data,x=T,y=T,surv=T)
data$p_lp <- predict(cox_GAP, data, type="lp")
c_harrell_GAP <- (cph(Surv(time,status)~p_lp, data=data,x=TRUE,y=TRUE)$stats["Dxy"]+1)/2
c_harrell_GAP


# internal validation
library(MASS)
N_bootstrap <- 1000      
c_harrell_resample <- 0      
c_harrell_original <- 0    
for (i in 1:N_bootstrap){
  data.train <- data[sample(1:nrow(data), replace=TRUE),]
  Outcome <- "Surv(time,status)"
  CandidateVariables <- c("age","RA","Immunosuppressive.agents","DLCO","RVD","RAA","AST","ALB","Honeycombing")
  Formula <- formula(paste(paste(Outcome,"~", collapse=" "), 
                           paste(CandidateVariables, collapse=" + ")))
  
  model.full <- coxph(Formula, data=data.train,x=TRUE)
  
  model.train <- stepAIC(model.full, direction="both")
  
  
  c_harrell_resample[i] <- (cph(Surv(time,status)~p_lp, data=data.train,x=TRUE,y=TRUE)$stats["Dxy"]+1)/2
  data.test <- data
  data.test$p_lp <- predict(model.train, data.test, type="lp")
  c_harrell_original[i] <- (cph(Surv(time,status)~p_lp, data=data.test,x=TRUE,y=TRUE)$stats["Dxy"]+1)/2
}

c_harrell_optimism <- mean(c_harrell_resample - c_harrell_original)
c_harrell_optimism
c_harrell_nomogram - c_harrell_optimism



# Calibration plots
nomogram_3y <- cph(Surv(time,status)~age+RA+Immunosuppressive.agents+DLCO+RVD+RAA+AST+ALB+Honeycombing,x=T, y=T, surv=T, time.inc = 36, data=data)
cal3_nomogram <- calibrate(nomogram_3y, cmethod="KM", method="boot", u=36, m= 168, B= 1000,conf.int="TRUE")
plot(cal3_nomogram,xlab="Predicted 3 Years Survival",ylab="Fraction Surviving 3 years")

nomogram_5y <- cph(Surv(time,status)~age+RA+Immunosuppressive.agents+DLCO+RVD+RAA+AST+ALB+Honeycombing,x=T, y=T, surv=T, time.inc = 60, data=data)
cal5_nomogram <- calibrate(nomogram_5y, cmethod="KM", method="boot", u=60, m= 168, B= 1000,conf.int="TRUE")
plot(cal5_nomogram,xlab="Predicted 5 Years Survival",ylab="Fraction Surviving 5 years")

GAP_3y <- cph(Surv(time,status)~GAP,x=T, y=T, surv=T, time.inc = 36, data=data)
cal3_GAP <- calibrate(GAP_3y, cmethod="KM", method="boot", u=36, m= 168, B= 1000,conf.int="TRUE")
plot(cal3_GAP,xlab="Predicted 3 Years Survival",ylab="Fraction Surviving 3 years")

GAP_5y <- cph(Surv(time,status)~GAP,x=T, y=T, surv=T, time.inc = 60, data=data)
cal5_GAP <- calibrate(GAP_5y, cmethod="KM", method="boot", u=60, m= 168, B= 1000,conf.int="TRUE")
plot(cal5_GAP,xlab="Predicted 5 Years Survival",ylab="Fraction Surviving 5 years")


# net reclassification improvement (NRI), integrated discrimination improvement (IDI) M1 is IDI, M2 is NRI
library(survival)
library(survC1)
library(survIDINRI)
# 3year
model.1 <- coxph(Surv(time,status)~age+RA+Immunosuppressive.agents+DLCO+RVD+RAA+AST+ALB+Honeycombing,data=data,x=TRUE)
m11 <- predict(model.1, data=data, type="lp")
model.2 <- coxph(Surv(time,status)~GAP,data=data,x=TRUE)
m22 <- predict(model.2, data=data, type="lp")
IDI<-IDI.INF(data[,c("time","status")],m22, m11, 36, npert = 300, npert.rand = NULL, seed1 = NULL, alpha = 0.05)
IDI.INF.OUT(IDI)
# 5year
IDI<-IDI.INF(data[,c("time","status")],m22, m11, 60, npert = 300, npert.rand = NULL, seed1 = NULL, alpha = 0.05)
IDI.INF.OUT(IDI)


# likelihood-ratio test
TN.GAP <- cph(Surv(time,status)~ GAP, 
                data=data, na.action=na.omit )
TNC.nomogram <- cph(Surv(time,status)~ age+RA+Immunosuppressive.agents+DLCO+RVD+RAA+AST+ALB+Honeycombing, 
                   data=data, na.action=na.omit )
TNC.combine <- cph(Surv(time,status)~ GAP+age+RA+Immunosuppressive.agents+DLCO+RVD+RAA+AST+ALB+Honeycombing, 
                 data=data, na.action=na.omit )
TN.GAP
TNC.nomogram
TNC.combine
TNC1 <- lrtest(TN.GAP, TNC.combine)
TNC1
TNC2 <- lrtest(TNC.nomogram, TNC.combine)
TNC2


# decision curve analysis (DCA)
library(rms)
setwd("D:/R/pac")
source("stdca.R")
# 3 year
coxmod <- cph(Surv(time,status)~age+RA+Immunosuppressive.agents+DLCO+RVD+RAA+AST+ALB+Honeycombing,x=T, y=T, surv=T, time.inc = 36,data=data)
data$our_model <- c(1 - (summary(survfit(coxmod,newdata=data), times=36)$surv))
coxmod1 <- cph(Surv(time,status)~GAP,x=T, y=T, surv=T, time.inc = 36,data=data)
data$GAP <- c(1 - (summary(survfit(coxmod1,newdata=data), times=36)$surv))
stdca(data=data, outcome="status", ttoutcome="time", timepoint=36, predictors=c("our_model","GAP"), xstop=0.5, smooth=TRUE)
# 5 year
coxmod <- cph(Surv(time,status)~age+RA+Immunosuppressive.agents+DLCO+RVD+RAA+AST+ALB+Honeycombing,x=T, y=T, surv=T, time.inc = 60,data=data)
data$our_model <- c(1 - (summary(survfit(coxmod,newdata=data), times=60)$surv))
coxmod1 <- cph(Surv(time,status)~GAP,x=T, y=T, surv=T, time.inc = 60,data=data)
data$GAP <- c(1 - (summary(survfit(coxmod1,newdata=data), times=60)$surv))
stdca(data=data, outcome="status", ttoutcome="time", timepoint=60, predictors=c("our_model","GAP"), xstop=0.5, smooth=TRUE)


# All-cause mortality among 675 CTD-ILD patients
library(Matrix)
library(glmnet)
library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)
A <- read_excel("./data/A.xlsx")
fit <- survfit(Surv(time,status)~1,data=A)
res.sum <- surv_summary((fit))
res.sum
ggsurvplot(fit,
           pval=T,conf.int = T,
           risk.table=T,
           risik.table.col="strata",
           linetype="strata",
           surv.median.line="hv",
           ggtheme=theme_bw(),
           palette=c("#E7B800","#2E9FDF"),
           fun = "event")
