#library(magrittr)



#--------------------------
# Manipulating data

#setwd("C:\Users\paulw\Documents\GitHub/EBPGrouped")
load("export_chiapas.RData")

survey_data_hh <- export_chiapas_surv
census_hh <- export_chiapas_cens

#------------------------------------------------------------------------------
#
# Data preperation
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
#census_hh<-census_hh[!(census_hh$inglabpc>100000),]
survey_data_hh<-survey_data_hh[!(survey_data_hh$inglabpc<1),]
#survey_data_hh<-survey_data_hh[!(survey_data_hh$inglabpc>100000),]
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

# Include only the reslevant variables
#census<-census_hh%>%select(inglabpc,bienes,actcom,pcocup,pcpering,jnived,clase_hog,cluster_id)
# Define data
census<-census_hh
survey_data<-survey_data_hh
rm("census_hh","survey_data_hh")

census_syn<-census
rm("census")

# 6) Daten fuer meine Syntax kompatibel machen
census_syn$y <- census_syn$inglabpc
census_syn$clusterid <- census_syn$cluster_id

survey_data$y <- survey_data$inglabpc
survey_data$clusterid <- survey_data$cluster_id


pop_data <- census_syn
samp_data <- survey_data

# 7) Das Modell testen
#library(MuMIn)

#Models
#mod <- lm(y~bienes+actcom+pcocup+pcpering+jnived+clase_hog,data=samp_data) 
#mod_mix <- lme(y~bienes+actcom+pcocup+pcpering+jnived+clase_hog, random=~1|cluster_id,method="REML",data = samp_data)
#mod_mix2 <- lme(y~bienes+actcom+pcocup+jnived+clase_hog, random=~1|clusterid,method="REML",data = samp_data)
#mod_mix2 <- lmer(y~bienes+actcom+pcocup+jnived+clase_hog + (1|clusterid),data = samp_data)


#Summary
#summary(mod)
#summary(mod_mix)
#summary(mod_mix2)

#R^2
#r.squaredGLMM(mod)
#r.squaredGLMM(mod_mix)
#r.squaredGLMM(mod_mix2)


#
# Box Cox Modell testen
#

#formel <- as.formula(y ~ bienes+actcom+pcocup+jnived+clase_hog)
#formula <- "y ~ bienes+actcom+pcocup+jnived+clase_hog+(1|clusterid)"

#source("BoxCox_function_final.R")
#sampBox <- boxcox(dat=samp_data, inverse = FALSE, formula = as.formula(y ~ bienes+actcom+pcocup+jnived+clase_hog))
#samp_dataBox <- samp_data
#samp_data$y <- sampBox[[1]]
#m <- sampBox[[2]]
#lambda <- sampBox[[3]]
#itEM <- NULL

#lmreg <- lmer(formula, data=samp_data)
#summary(lmreg)
#r.squaredGLMM(lmreg)

#ICC
#icc <- function(u,e){u/(u+e)}

#se <- mod_mix$sigma^2
#su <- as.numeric(VarCorr(mod_mix)[1,1])
#icc(u=su,e=se)

#classes <- c(1,50,100,200,400,600,1000,1500,2000,3000,4000,5500,8000,12000,Inf)
#samp_data$yclass <- cut(samp_data$y, classes) 
#table(samp_data$yclass)
#formula <- "yclass ~ bienes+actcom+pcocup+jnived+clase_hog+(1|clusterid)"
#library(smicd)
#sem <- semLme(formula = formula, classes = classes, burnin = 10, samples = 40, 
 #             data = samp_data, trafo = "bc")
#summary(sem)


rm("census_syn","export_chiapas_cens", "export_chiapas_surv", "survey_data")

#------------------------------------------------------------------------------
#
# EBP emdi without area 100 and 101
#"y ~ bienes+actcom+pcocup+jnived+clase_hog+(1|clusterid)"
#------------------------------------------------------------------------------
#samp_data <- samp_data[samp_data$clusterid!=101,]
#samp_data <- samp_data[samp_data$clusterid!=100,]

#pop_data <- pop_data[pop_data$clusterid!=101,]
#pop_data <- pop_data[pop_data$clusterid!=100,]

#library(emdi)

#direct_estim <- direct("y", smp_data = samp_data, smp_domains = "clusterid")

#ebp_estim <- ebp(y ~ bienes + actcom + pcocup + jnived + clase_hog,
                 #pop_data = pop_data, pop_domains = "clusterid", 
                 #smp_data = samp_data, smp_domains = "clusterid", transformation = "box.cox")

#compare_plot(direct_estim, ebp_estim, indicator = "Mean")

#summary(ebp_estim)


#------------------------------------------------------------------------------
#
# EBP
#
#------------------------------------------------------------------------------
library(ineq)
#library(intReg, lib.loc="H:/Doktor2/libWin")
library(pps)
#library(MASS, lib.loc="H:/Doktor2/libWin")
library(nlme)
#library(sp, lib.loc="H:/Doktor2/libWin")
#library(spgwr, lib.loc="H:/Doktor2/libWin")
#library(spdep, lib.loc="H:/Doktor2/libWin")
#library(mnormt, lib.loc="H:/Doktor2/libWin")
#library("actuar", lib.loc="/home/paulwalter/R/x86_64-pc-linux-gnu-library/3.2/")
library(lme4)
#library(FNN, lib.loc="H:/Doktor2/libWin")
library(truncnorm)
#library(psych, lib.loc="H:/Doktor2/libWin")
#library(readstata13, lib.loc="H:/Doktor2/libWin")
#library(foreign, lib.loc="H:/Doktor2/libWin")

#library(lqmm, lib.loc="H:/Doktor2/libWin")
#library(reshape, lib.loc="H:/Doktor2/libWin")
#library(dplyr, lib.loc="H:/Doktor2/libWin")
#library(laeken, lib.loc="H:/Doktor2/libWin")
#library(quantreg, lib.loc="H:/Doktor2/libWin")
library(operator.tools)
#library(operator.tools)
library(formula.tools)
#library(formula.tools)
#library(mitml, lib.loc="H:/Doktor2/libWin")
library(MuMIn)

#classes <- c(1,50,100,200,400,600,1000,1500,2250,3000,4000,5500,8000,12000,Inf)


## Hier die Klassierung ändern!
classes <- c(1,50,100,200,400,600,1000,1500,2000,3000,4000,5500,8000,12000,Inf)
# classes <- c(1,100,400,1000,2000,4000,8000,12000,Inf)
#classes <- c(1,400,2000,8000,Inf)

samp_data$yclass <- cut(samp_data$y, classes) 
plot(samp_data$yclass)
table(samp_data$yclass)

# Zu setzende Parameter


#burnin <- c(4)  #40
#samples <- c(20) #200 
#loop <- 10 #100
#boot <- 10 #100 
#terminate <- 0.005
#MSE <- TRUE
#method = c("LME_box", "SEM_box")

formel <- as.formula(y ~ bienes+actcom+pcocup+jnived+clase_hog)
formula <- "y ~ bienes+actcom+pcocup+jnived+clase_hog+(1|clusterid)"

pop_data$idD <- pop_data$cluster_id
samp_data$idD <- samp_data$cluster_id


source("MSE_function_final_rueck.R")

# LME Box-Cox
results_lme_box <- mse(formulaFix=formel, formula=formula, population=pop_data , 
               sample=samp_data, L=5, threshold=0.5, MSE=FALSE , B=5,
               burnin=5, samples=20, classes, method="LME", terminate = 0.01, transformation ="box"
        )

# SEM Box-Cox
results_sem_box <- mse(formulaFix=formel, formula=formula, population=pop_data , 
               sample=samp_data, L=5, threshold=0.5, MSE=FALSE , B=5,
               burnin=5, samples=20, classes, method="SEM", terminate = 0.01, transformation ="box"
)


#--------------------------
# Some results LME Box-Cox
mean <- as.vector(results_lme_box[[1]][[1]]$mean)
hcr <- as.vector(results_lme_box[[1]][[1]]$hcr)
pgap <- as.vector(results_lme_box[[1]][[1]]$pgap)

# Some results SEM Box-Cox
mean <- as.vector(results_sem_box[[1]][[1]]$mean)
hcr <- as.vector(results_sem_box[[1]][[1]]$hcr)
pgap <- as.vector(results_sem_box[[1]][[1]]$pgap)

