#Author: Paul Walter
#Mail: paul.walter@fu-berlin.de


library(ineq)
library(pps)
library(nlme)
library(lme4)
library(truncnorm)
library(operator.tools)
library(formula.tools)
library(MuMIn)

#--------------------------
# Import data

load("data/export_chiapas.RData")
survey_data_hh <- export_chiapas_surv
census_hh <- export_chiapas_cens

#------------------------------------------------------------------------------
#
# Data prep
#------------------------------------------------------------------------------

# 1) Create AREA IDs for Sample and Census
census_hh<-data.frame(census_hh,cluster_id=as.numeric(census_hh$mun))
census_hh<-data.frame(census_hh, rururb2=as.numeric(census_hh$rururb=="Urbano"))
survey_data_hh<-data.frame(survey_data_hh,cluster_id=as.numeric(survey_data_hh$mun))


# 2) Sorting of the data
census_hh<-census_hh[order(census_hh$cluster_id), ] 
survey_data_hh<-survey_data_hh[order(survey_data_hh$cluster_id), ] 

# 4) Remove 0 income or below
census_hh<-census_hh[!(census_hh$inglabpc<1),]
survey_data_hh<-survey_data_hh[!(survey_data_hh$inglabpc<1),]
summary(census_hh$inglabpc)
summary(survey_data_hh$inglabpc)

# 5) Remove NA for working model
census_hh <- census_hh[!(is.na(census_hh$inglabpc)),]
census_hh <- census_hh[!(is.na(census_hh$bienes)),]
census_hh <- census_hh[!(is.na(census_hh$actcom)),]
census_hh <- census_hh[!(is.na(census_hh$pcocup)),]
census_hh <- census_hh[!(is.na(census_hh$pcpering)),]
census_hh <- census_hh[!(is.na(census_hh$jnived)),]
census_hh <- census_hh[!(is.na(census_hh$clase_hog)),]

survey_data_hh <- survey_data_hh[!(is.na(survey_data_hh$inglabpc)),]
survey_data_hh <- survey_data_hh[!(is.na(survey_data_hh$bienes)),]
survey_data_hh <- survey_data_hh[!(is.na(survey_data_hh$actcom)),]
survey_data_hh <- survey_data_hh[!(is.na(survey_data_hh$pcocup)),]
survey_data_hh <- survey_data_hh[!(is.na(survey_data_hh$pcpering)),]
survey_data_hh <- survey_data_hh[!(is.na(survey_data_hh$jnived)),]
survey_data_hh <- survey_data_hh[!(is.na(survey_data_hh$clase_hog)),]

census<-census_hh
survey_data<-survey_data_hh
rm("census_hh","survey_data_hh")

census_syn<-census
rm("census")

# 6) Rename data to make it compatible to the custom syntax
census_syn$y <- census_syn$inglabpc
census_syn$clusterid <- census_syn$cluster_id

survey_data$y <- survey_data$inglabpc
survey_data$clusterid <- survey_data$cluster_id

pop_data <- census_syn
samp_data <- survey_data


rm("census_syn","export_chiapas_cens", "export_chiapas_surv", "survey_data")



## Group data (if not grouped by the data provider already)
classes <- c(1,50,100,200,400,600,1000,1500,2000,3000,4000,5500,8000,12000,Inf)
samp_data$yclass <- cut(samp_data$y, classes) 
plot(samp_data$yclass)
table(samp_data$yclass)

# Formula
formel <- as.formula(y ~ bienes+actcom+pcocup+jnived+clase_hog)
formula <- "y ~ bienes+actcom+pcocup+jnived+clase_hog+(1|clusterid)"

# Rename of domain ID
pop_data$idD <- pop_data$cluster_id
samp_data$idD <- samp_data$cluster_id


source("helper/estimation_main.R")

# LME Box-Cox
results_lme_box <- EBP(formulaFix=formel, 
                       formula=formula, 
                       population=pop_data , 
                       sample=samp_data, 
                       L=5,  # for EBP
                       threshold=0.5, #for poverty estimators  
                       MSE=TRUE , 
                       B=5, # for MSE
                       burnin=5, # for SEM
                       samples=20,  # for SEM
                       classes, 
                       method="LME", #LME or SEM
                       transformation ="box" #box, log or none
        )

# SEM Box-Cox
results_sem_box <- EBP(formulaFix=formel, formula=formula, population=pop_data , 
               sample=samp_data, L=5, threshold=0.5, MSE=TRUE , B=5,
               burnin=5, samples=20, classes, method="SEM",  transformation ="box"
)


#--------------------------
# Some results LME Box-Cox
mean_lme <- as.vector(results_lme_box[[1]][[1]]$mean)
mean_lme_mse <- results_lme_box$MSE[,"mean"]
hcr_lme <- as.vector(results_lme_box[[1]][[1]]$hcr)
hcr_lme_mse <- results_lme_box$MSE[,"hcr"]
pgap_lme <- as.vector(results_lme_box[[1]][[1]]$pgap)
pgap_lme_mse <- results_lme_box$MSE[,"pgap"]

# Some results SEM Box-Cox
mean_sem <- as.vector(results_sem_box[[1]][[1]]$mean)
mean_sem_mse <- results_lme_box$MSE[,"mean"]
hcr_sem <- as.vector(results_sem_box[[1]][[1]]$hcr)
hcr_sem_mse <- results_sem_box$MSE[,"hcr"]
pgap_sem <- as.vector(results_sem_box[[1]][[1]]$pgap)
pgap_sem_mse <- results_sem_box$MSE[,"hcr"]

