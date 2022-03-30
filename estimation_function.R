# Simulation

sim <- function(burnin_initial, samples_initial ,  l, b,t, MSE=FALSE, method,
                #method= c("SEM", "SEM_log", "SEM_box", "EM", "EM_log", "EM_box", "LME", "LME_log", "LME_box"),  
                samp_data, pop_data, formel, formula){
  
  # Server
  #load("Data.RData")
  
  
  # Performance measures
#  funIndicators <- function(x , threshold) {
#    c(mean = mean(x), 
#      gini = ineq(x),
#      hcr = mean(x<threshold),
#      quant10 = unname(quantile(x, probs = c(0.1))),
#      quant25 = unname(quantile(x, probs = c(0.25))),
#      quant50 = unname(quantile(x, probs = c(0.50))),
#      quant75 = unname(quantile(x, probs = c(0.75))),
#      quant90 = unname(quantile(x, probs = c(0.90))),
#      pgap = mean((x<threshold)*(threshold-x)/threshold),
#      qsr = sum(x[(x>quantile(x,0.8))])/sum(x[(x<quantile(x,0.2))]))
 # }
  
  
  #area_SEM_estimates <- cbind(EBP[[1]]$mean, EBP[[1]]$gini, EBP[[1]]$hcr, EBP[[1]]$quant10, EBP[[1]]$quant25, 
  #                            EBP[[1]]$quant75, EBP[[1]]$quant90)
  
  pop_estimates <- vector("list", length = 1)
  area_pop_estimates <- vector("list", length = 1) 
  
  #for (i in 1:length(Pop) ){
  threshold <- 0.6*median(pop_data$y)
  pop_estimates <- funIndicators(pop_data$y, threshold = threshold)
  pop_data$clusterid <- pop_data$idD
  samp_data$clusterid <- samp_data$idD
  #}
  
  #for (i in 1:length(Pop)) {
  area_pop_estimates <- by(pop_data, pop_data$clusterid, function(x){funIndicators(x$y, threshold = threshold)})
  #}
  
  
  #------------------------------------------------------------------------------
  
  # Speichern NEU
  
  #------------------------------------------------------------------------------
  #list.names <- c(method, "True")
  
  
  list.names <- c("LME", "LME_log", "LME_box", "MID", "SEM", "SEM_log", "SEM_box", "EM", "EM_log", "EM_box" , "True")
  list.results <- vector("list", length(list.names))
  names(list.results) <- list.names
  
  list.names <- c("true.mean", "true.gini", "true.hcr", "true.quant10", "true.quant25","true.quant50", "true.quant75", "true.quant90", "true.pgap", "true.qsr")
  list.True <- vector("list", length = length(pop_estimates))
  names(list.True) <- list.names
  help <- matrix(0, nrow = 1, ncol = length(unique(pop_data$clusterid)))
  
  for (i in 1:length(list.True)) {
    list.True[[i]] <- help
  }
  
  list.results$True <- list.True
  
  list.names <- c("est.mean", "est.gini", "est.hcr", "est.quant10", "est.quant25","est.quant50", "est.quant75", "est.quant90", "est.pgap", "est.qsr",
                  "MSE.mean", "MSE.gini", "MSE.hcr", "MSE.quant10", "MSE.quant25","MSE.quant50","MSE.quant75", "MSE.quant90", "MSE.pgap", "MSE.qsr")
  list.indicators <- vector("list", length = length(pop_estimates)*2)
  names(list.indicators)<- list.names
  help <- matrix(0, nrow = 1, ncol = length(unique(pop_data$clusterid)))
  
  for (i in 1:length(list.indicators)) {
    list.indicators[[i]] <- help
  }
  
  for (i in 1:(length(list.results)-1)) {
    list.results[[i]]<-list.indicators
  }
  
  for (i in 1:1) {
    for (j in 1:length(unique(pop_data$clusterid))){
      list.results$True$true.mean[i,j] <- area_pop_estimates[[j]]["mean"]
      list.results$True$true.gini[i,j] <- area_pop_estimates[[j]]["gini"]
      list.results$True$true.hcr[i,j] <- area_pop_estimates[[j]]["hcr"]
      list.results$True$true.quant10[i,j] <- area_pop_estimates[[j]]["quant10"]
      list.results$True$true.quant25[i,j] <- area_pop_estimates[[j]]["quant25"]
      list.results$True$true.quant50[i,j] <- area_pop_estimates[[j]]["quant50"]
      list.results$True$true.quant75[i,j] <- area_pop_estimates[[j]]["quant75"]
      list.results$True$true.quant90[i,j] <- area_pop_estimates[[j]]["quant90"]
      list.results$True$true.pgap[i,j] <- area_pop_estimates[[j]]["pgap"]
      list.results$True$true.qsr[i,j] <- area_pop_estimates[[j]]["qsr"]
    }
  }
  
  
  #==============================================================================
  # Save Lambda from each Simulation run
  
  
  LME_Lambda_Sim <- vector(mode="numeric", length = 1)
  EM_Lambda_Sim <- vector(mode="numeric", length = 1)
  SEM_Lambda_Sim <- vector(mode="numeric", length = 1)
  
  #==============================================================================
  # Save predicted y population for one particular simulation run
  
  SEM_EBP_pop <- vector()
  SEM_EBP_pop_log <- vector()
  SEM_EBP_pop_box <- vector()
  EM_EBP_pop <- vector()
  EM_EBP_pop_log <- vector()
  EM_EBP_pop_box <- vector()
  LME_EBP_pop <- vector()
  LME_EBP_pop_log <- vector()
  LME_EBP_pop_box <- vector()
  MID_EBP_pop <- vector()
  
  #==============================================================================
  # Save iterations EM 
  itEM <- vector()
  itEMLog <- vector()
  itEMBox <- vector()
  
  true_y <- pop_data$y
  
  for (i in 1:1){
    
    #------------------------------------------------------------------------------
    
    # Drawing of samples
    #------------------------------------------------------------------------------
    
    #if (exists("balanced")==TRUE) {
    #if (balanced==TRUE) {s<-stratsrs(pop_data$clusterid,rep(SampleSetup$snunits,DataSetup$nclusters))}
    #if (balanced==FALSE) {s<-stratsrs(pop_data$clusterid,sample_size)}
    #}
    
    #if (exists("balanced")==FALSE) { s<-stratsrs(pop_data$clusterid,sample_size)}
    
    #samp_data <- Pop[[no]]
    #pop_data <- attr(Pop[[no]], "pop")
    
    # Initialize function calls EM, SEM
    #formel <- as.formula(y~x)
    #formula <- "y ~ x+(1|clusterid)"
    burnin <- burnin_initial
    samples <- samples_initial
    
    #==============================================================================
    # Apply function MID
    
    if ("MID" %in% method) {
      source("MSE_function_final_rueck.R")
      EBP <- mse(formulaFix=formel,  formula=formula, population=pop_data , sample=samp_data, L=l, threshold=threshold, 
                 MSE=MSE, B=b, burnin=burnin, samples=samples, classes=classes , method="MID", transformation = "none")
      
      MID_EBP_pop <-EBP[[1]]$y  
      
      # Neue art zu speichern
      list.results$MID$est.mean[i,] <- EBP[[1]][[1]]$mean
      list.results$MID$est.gini[i,] <- EBP[[1]][[1]]$gini
      list.results$MID$est.hcr [i,] <- EBP[[1]][[1]]$hcr 
      list.results$MID$est.quant10[i,] <- EBP[[1]][[1]]$quant10
      list.results$MID$est.quant25[i,] <- EBP[[1]][[1]]$quant25
      list.results$MID$est.quant50[i,] <- EBP[[1]][[1]]$quant50
      list.results$MID$est.quant75[i,] <- EBP[[1]][[1]]$quant75
      list.results$MID$est.quant90[i,] <- EBP[[1]][[1]]$quant90
      list.results$MID$est.pgap[i,] <- EBP[[1]][[1]]$pgap
      list.results$MID$est.qsr[i,] <- EBP[[1]][[1]]$qsr
      
      if (MSE==TRUE){
        list.results$MID$MSE.mean[i,] <- EBP$MSE[,"mean"]
        list.results$MID$MSE.gini[i,] <- EBP$MSE[,"gini"]
        list.results$MID$MSE.hcr[i,] <- EBP$MSE[,"hcr"]
        list.results$MID$MSE.quant10[i,] <- EBP$MSE[,"quant10"]
        list.results$MID$MSE.quant25[i,] <- EBP$MSE[,"quant25"]
        list.results$MID$MSE.quant50[i,] <- EBP$MSE[,"quant50"]
        list.results$MID$MSE.quant75[i,] <- EBP$MSE[,"quant75"]
        list.results$MID$MSE.quant90[i,] <- EBP$MSE[,"quant90"]
        list.results$MID$MSE.pgap[i,] <- EBP$MSE[,"pgap"]
        list.results$MID$MSE.qsr[i,] <- EBP$MSE[,"qsr"]
      }
    }
    
    #==============================================================================
    # Apply function SEM
    
    # EBP estimates (neuer Aufruf mit MSE)
    if ("SEM" %in% method) {
      source("MSE_function_final_rueck.R")
      EBP <- mse(formulaFix=formel,  formula=formula, population=pop_data , sample=samp_data, L=l, threshold=threshold, 
                 MSE=MSE, B=b, burnin=burnin, samples=samples, classes=classes , method="SEM", transformation = "none")
      
      SEM_EBP_pop <-EBP[[1]]$y  
      
      # Neue art zu speichern
      list.results$SEM$est.mean[i,] <- EBP[[1]][[1]]$mean
      list.results$SEM$est.gini[i,] <- EBP[[1]][[1]]$gini
      list.results$SEM$est.hcr [i,] <- EBP[[1]][[1]]$hcr 
      list.results$SEM$est.quant10[i,] <- EBP[[1]][[1]]$quant10
      list.results$SEM$est.quant25[i,] <- EBP[[1]][[1]]$quant25
      list.results$SEM$est.quant50[i,] <- EBP[[1]][[1]]$quant50
      list.results$SEM$est.quant75[i,] <- EBP[[1]][[1]]$quant75
      list.results$SEM$est.quant90[i,] <- EBP[[1]][[1]]$quant90
      list.results$SEM$est.pgap[i,] <- EBP[[1]][[1]]$pgap
      list.results$SEM$est.qsr[i,] <- EBP[[1]][[1]]$qsr
      
      if (MSE==TRUE){
        list.results$SEM$MSE.mean[i,] <- EBP$MSE[,"mean"]
        list.results$SEM$MSE.gini[i,] <- EBP$MSE[,"gini"]
        list.results$SEM$MSE.hcr[i,] <- EBP$MSE[,"hcr"]
        list.results$SEM$MSE.quant10[i,] <- EBP$MSE[,"quant10"]
        list.results$SEM$MSE.quant25[i,] <- EBP$MSE[,"quant25"]
        list.results$SEM$MSE.quant50[i,] <- EBP$MSE[,"quant50"]
        list.results$SEM$MSE.quant75[i,] <- EBP$MSE[,"quant75"]
        list.results$SEM$MSE.quant90[i,] <- EBP$MSE[,"quant90"]
        list.results$SEM$MSE.pgap[i,] <- EBP$MSE[,"pgap"]
        list.results$SEM$MSE.qsr[i,] <- EBP$MSE[,"qsr"]
      }
    }
    
    #==============================================================================
    # Apply function SEM Log
    
    # EBP estimates (neuer Aufruf mit MSE)
    if ("SEM_log" %in% method) {
      source("MSE_function_final_rueck.R")
      EBP <- mse(formulaFix=formel,  formula=formula, population=pop_data , sample=samp_data, L=l, threshold=threshold, 
                 MSE=MSE, B=b, burnin=burnin, samples=samples, classes=classes , method="SEM" , transformation="log")
      
      SEM_EBP_pop_log <-EBP[[1]]$y  
      
      list.results$SEM_log$est.mean[i,] <- EBP[[1]][[1]]$mean
      list.results$SEM_log$est.gini[i,] <- EBP[[1]][[1]]$gini
      list.results$SEM_log$est.hcr [i,] <- EBP[[1]][[1]]$hcr 
      list.results$SEM_log$est.quant10[i,] <- EBP[[1]][[1]]$quant10
      list.results$SEM_log$est.quant25[i,] <- EBP[[1]][[1]]$quant25
      list.results$SEM_log$est.quant50[i,] <- EBP[[1]][[1]]$quant50
      list.results$SEM_log$est.quant75[i,] <- EBP[[1]][[1]]$quant75
      list.results$SEM_log$est.quant90[i,] <- EBP[[1]][[1]]$quant90
      list.results$SEM_log$est.pgap[i,] <- EBP[[1]][[1]]$pgap
      list.results$SEM_log$est.qsr[i,] <- EBP[[1]][[1]]$qsr
      
      if (MSE==TRUE){
        list.results$SEM_log$MSE.mean[i,] <- EBP$MSE[,"mean"]
        list.results$SEM_log$MSE.gini[i,] <- EBP$MSE[,"gini"]
        list.results$SEM_log$MSE.hcr[i,] <- EBP$MSE[,"hcr"]
        list.results$SEM_log$MSE.quant10[i,] <- EBP$MSE[,"quant10"]
        list.results$SEM_log$MSE.quant25[i,] <- EBP$MSE[,"quant25"]
        list.results$SEM_log$MSE.quant50[i,] <- EBP$MSE[,"quant50"]
        list.results$SEM_log$MSE.quant75[i,] <- EBP$MSE[,"quant75"]
        list.results$SEM_log$MSE.quant90[i,] <- EBP$MSE[,"quant90"]
        list.results$SEM_log$MSE.pgap[i,] <- EBP$MSE[,"pgap"]
        list.results$SEM_log$MSE.qsr[i,] <- EBP$MSE[,"qsr"]
      }
    }
    
    #==============================================================================
    # Apply function SEM BoxCox
    
    
    # EBP estimates (neuer Aufruf mit MSE)
    if ("SEM_box" %in% method){
      source("MSE_function_final_rueck.R")
      EBP <- mse(formulaFix=formel,  formula=formula, population=pop_data , sample=samp_data, L=l, threshold=threshold, 
                 MSE=MSE, B=b, burnin=burnin, samples=samples, classes=classes , method="SEM" , transformation="box")
      
      SEM_EBP_pop_box <-EBP[[1]]$y  
      
      list.results$SEM_box$est.mean[i,] <- EBP[[1]][[1]]$mean
      list.results$SEM_box$est.gini[i,] <- EBP[[1]][[1]]$gini
      list.results$SEM_box$est.hcr [i,] <- EBP[[1]][[1]]$hcr 
      list.results$SEM_box$est.quant10[i,] <- EBP[[1]][[1]]$quant10
      list.results$SEM_box$est.quant25[i,] <- EBP[[1]][[1]]$quant25
      list.results$SEM_box$est.quant50[i,] <- EBP[[1]][[1]]$quant50
      list.results$SEM_box$est.quant75[i,] <- EBP[[1]][[1]]$quant75
      list.results$SEM_box$est.quant90[i,] <- EBP[[1]][[1]]$quant90
      list.results$SEM_box$est.pgap[i,] <- EBP[[1]][[1]]$pgap
      list.results$SEM_box$est.qsr[i,] <- EBP[[1]][[1]]$qsr
      if (MSE==FALSE){SEM_Lambda_Sim[i] <- EBP$lambda}
      
      if (MSE==TRUE){
        list.results$SEM_box$MSE.mean[i,] <- EBP$MSE[,"mean"]
        list.results$SEM_box$MSE.gini[i,] <- EBP$MSE[,"gini"]
        list.results$SEM_box$MSE.hcr[i,] <- EBP$MSE[,"hcr"]
        list.results$SEM_box$MSE.quant10[i,] <- EBP$MSE[,"quant10"]
        list.results$SEM_box$MSE.quant25[i,] <- EBP$MSE[,"quant25"]
        list.results$SEM_box$MSE.quant50[i,] <- EBP$MSE[,"quant50"]
        list.results$SEM_box$MSE.quant75[i,] <- EBP$MSE[,"quant75"]
        list.results$SEM_box$MSE.quant90[i,] <- EBP$MSE[,"quant90"]
        list.results$SEM_box$MSE.pgap[i,] <- EBP$MSE[,"pgap"]
        list.results$SEM_box$MSE.qsr[i,] <- EBP$MSE[,"qsr"]
        SEM_Lambda_Sim[i]<- EBP[[3]]$lambda
      }
    }
    #==============================================================================
    # Apply function EM
    
    # EBP estimates (neuer Aufruf mit MSE)
    if ("EM" %in% method) {
      source("MSE_function_final_rueck.R")
      EBP <- mse(formulaFix=formel,  formula=formula, population=pop_data , sample=samp_data, L=l, threshold=threshold, 
                 MSE=MSE ,B=b, burnin=burnin, samples=samples, classes=classes , method="EM" , terminate = t, transformation = "none")
      
      EM_EBP_pop <-EBP[[1]]$y  
      
      
      list.results$EM$est.mean[i,] <- EBP[[1]][[1]]$mean
      list.results$EM$est.gini[i,] <- EBP[[1]][[1]]$gini
      list.results$EM$est.hcr [i,] <- EBP[[1]][[1]]$hcr 
      list.results$EM$est.quant10[i,] <- EBP[[1]][[1]]$quant10
      list.results$EM$est.quant25[i,] <- EBP[[1]][[1]]$quant25
      list.results$EM$est.quant50[i,] <- EBP[[1]][[1]]$quant50
      list.results$EM$est.quant75[i,] <- EBP[[1]][[1]]$quant75
      list.results$EM$est.quant90[i,] <- EBP[[1]][[1]]$quant90
      list.results$EM$est.pgap[i,] <- EBP[[1]][[1]]$pgap
      list.results$EM$est.qsr[i,] <- EBP[[1]][[1]]$qsr
      if (MSE==FALSE){itEM <- EBP$it}
      
      if (MSE==TRUE){
        list.results$EM$MSE.mean[i,] <- EBP$MSE[,"mean"]
        list.results$EM$MSE.gini[i,] <- EBP$MSE[,"gini"]
        list.results$EM$MSE.hcr[i,] <- EBP$MSE[,"hcr"]
        list.results$EM$MSE.quant10[i,] <- EBP$MSE[,"quant10"]
        list.results$EM$MSE.quant25[i,] <- EBP$MSE[,"quant25"]
        list.results$EM$MSE.quant50[i,] <- EBP$MSE[,"quant50"]
        list.results$EM$MSE.quant75[i,] <- EBP$MSE[,"quant75"]
        list.results$EM$MSE.quant90[i,] <- EBP$MSE[,"quant90"]
        list.results$EM$MSE.pgap[i,] <- EBP$MSE[,"pgap"]
        list.results$EM$MSE.qsr[i,] <- EBP$MSE[,"qsr"]
        itEM <- EBP[[3]]$it
      }
    }
    
    #==============================================================================
    # Apply function EM log
    
    # EBP estimates (neuer Aufruf mit MSE)
    if ("EM_log" %in% method){
      source("MSE_function_final_rueck.R")
      EBP <- mse(formulaFix=formel,  formula=formula, population=pop_data , sample=samp_data, L=l, threshold=threshold, MSE=MSE 
                 ,B=b, burnin=burnin, samples=samples, classes=classes , method="EM" , terminate = t, transformation = "log")
      
      EM_EBP_pop_log <-EBP[[1]]$y  
      
      list.results$EM_log$est.mean[i,] <- EBP[[1]][[1]]$mean
      list.results$EM_log$est.gini[i,] <- EBP[[1]][[1]]$gini
      list.results$EM_log$est.hcr [i,] <- EBP[[1]][[1]]$hcr 
      list.results$EM_log$est.quant10[i,] <- EBP[[1]][[1]]$quant10
      list.results$EM_log$est.quant25[i,] <- EBP[[1]][[1]]$quant25
      list.results$EM_log$est.quant50[i,] <- EBP[[1]][[1]]$quant50
      list.results$EM_log$est.quant75[i,] <- EBP[[1]][[1]]$quant75
      list.results$EM_log$est.quant90[i,] <- EBP[[1]][[1]]$quant90
      list.results$EM_log$est.pgap[i,] <- EBP[[1]][[1]]$pgap
      list.results$EM_log$est.qsr[i,] <- EBP[[1]][[1]]$qsr
      if (MSE==FALSE){itEMLog <- EBP$it}
      
      if (MSE==TRUE){
        list.results$EM_log$MSE.mean[i,] <- EBP$MSE[,"mean"]
        list.results$EM_log$MSE.gini[i,] <- EBP$MSE[,"gini"]
        list.results$EM_log$MSE.hcr[i,] <- EBP$MSE[,"hcr"]
        list.results$EM_log$MSE.quant10[i,] <- EBP$MSE[,"quant10"]
        list.results$EM_log$MSE.quant25[i,] <- EBP$MSE[,"quant25"]
        list.results$EM_log$MSE.quant50[i,] <- EBP$MSE[,"quant50"]
        list.results$EM_log$MSE.quant75[i,] <- EBP$MSE[,"quant75"]
        list.results$EM_log$MSE.quant90[i,] <- EBP$MSE[,"quant90"]
        list.results$EM_log$MSE.pgap[i,] <- EBP$MSE[,"pgap"]
        list.results$EM_log$MSE.qsr[i,] <- EBP$MSE[,"qsr"]
        itEMLog <- EBP[[3]]$it
      }
    }
    
    #==============================================================================
    # Apply function EM Box
    
    if ("EM_box" %in% method){
      source("MSE_function_final_rueck.R")
      EBP <- mse(formulaFix=formel,  formula=formula, population=pop_data , sample=samp_data, L=l, threshold=threshold, MSE=MSE 
                 ,B=b, burnin=burnin, samples=samples, classes=classes , method="EM" , terminate = t, transformation = "box")
      
      EM_EBP_pop_box <-EBP[[1]]$y  
      
      list.results$EM_box$est.mean[i,] <- EBP[[1]][[1]]$mean
      list.results$EM_box$est.gini[i,] <- EBP[[1]][[1]]$gini
      list.results$EM_box$est.hcr [i,] <- EBP[[1]][[1]]$hcr 
      list.results$EM_box$est.quant10[i,] <- EBP[[1]][[1]]$quant10
      list.results$EM_box$est.quant25[i,] <- EBP[[1]][[1]]$quant25
      list.results$EM_box$est.quant50[i,] <- EBP[[1]][[1]]$quant50
      list.results$EM_box$est.quant75[i,] <- EBP[[1]][[1]]$quant75
      list.results$EM_box$est.quant90[i,] <- EBP[[1]][[1]]$quant90
      list.results$EM_box$est.pgap[i,] <- EBP[[1]][[1]]$pgap
      list.results$EM_box$est.qsr[i,] <- EBP[[1]][[1]]$qsr
      if (MSE==FALSE) {EM_Lambda_Sim[i] <- EBP$lambda
      itEMBox <- EBP$it}
      
      if (MSE==TRUE){
        list.results$EM_box$MSE.mean[i,] <- EBP$MSE[,"mean"]
        list.results$EM_box$MSE.gini[i,] <- EBP$MSE[,"gini"]
        list.results$EM_box$MSE.hcr[i,] <- EBP$MSE[,"hcr"]
        list.results$EM_box$MSE.quant10[i,] <- EBP$MSE[,"quant10"]
        list.results$EM_box$MSE.quant25[i,] <- EBP$MSE[,"quant25"]
        list.results$EM_box$MSE.quant50[i,] <- EBP$MSE[,"quant50"]
        list.results$EM_box$MSE.quant75[i,] <- EBP$MSE[,"quant75"]
        list.results$EM_box$MSE.quant90[i,] <- EBP$MSE[,"quant90"]
        list.results$EM_box$MSE.pgap[i,] <- EBP$MSE[,"pgap"]
        list.results$EM_box$MSE.qsr[i,] <- EBP$MSE[,"qsr"]
        
        EM_Lambda_Sim[i] <- EBP[[3]]$lambda
        itEMBox <- EBP[[3]]$it
      }
    }
    
    #==============================================================================
    # Apply lm
    
    if ("LME" %in% method){
      source("MSE_function_final_rueck.R")
      EBP <- mse(formulaFix=formel,  formula=formula, population=pop_data , sample=samp_data, L=l, threshold=threshold, 
                 MSE=MSE, B=b , method="LME", transformation = "none")
      
      LME_EBP_pop <-EBP[[1]]$y 
      
      list.results$LME$est.mean[i,] <- EBP[[1]][[1]]$mean
      list.results$LME$est.gini[i,] <- EBP[[1]][[1]]$gini
      list.results$LME$est.hcr [i,] <- EBP[[1]][[1]]$hcr 
      list.results$LME$est.quant10[i,] <- EBP[[1]][[1]]$quant10
      list.results$LME$est.quant25[i,] <- EBP[[1]][[1]]$quant25
      list.results$LME$est.quant50[i,] <- EBP[[1]][[1]]$quant50
      list.results$LME$est.quant75[i,] <- EBP[[1]][[1]]$quant75
      list.results$LME$est.quant90[i,] <- EBP[[1]][[1]]$quant90
      list.results$LME$est.pgap[i,] <- EBP[[1]][[1]]$pgap
      list.results$LME$est.qsr[i,] <- EBP[[1]][[1]]$qsr
      
      if (MSE==TRUE){
        list.results$LME$MSE.mean[i,] <- EBP$MSE[,"mean"]
        list.results$LME$MSE.gini[i,] <- EBP$MSE[,"gini"]
        list.results$LME$MSE.hcr[i,] <- EBP$MSE[,"hcr"]
        list.results$LME$MSE.quant10[i,] <- EBP$MSE[,"quant10"]
        list.results$LME$MSE.quant25[i,] <- EBP$MSE[,"quant25"]
        list.results$LME$MSE.quant50[i,] <- EBP$MSE[,"quant50"]
        list.results$LME$MSE.quant75[i,] <- EBP$MSE[,"quant75"]
        list.results$LME$MSE.quant90[i,] <- EBP$MSE[,"quant90"]
        list.results$LME$MSE.pgap[i,] <- EBP$MSE[,"pgap"]
        list.results$LME$MSE.qsr[i,] <- EBP$MSE[,"qsr"]
      }
    }
    
    #==============================================================================
    # Apply lm log
    
    if ("LME_log" %in% method) {
      source("MSE_function_final_rueck.R")
      EBP <- mse(formulaFix=formel,  formula=formula, population=pop_data, sample=samp_data, L=l, threshold=threshold, 
                 MSE=MSE, B=b, method="LME", transformation="log")
      
      LME_EBP_pop_log <-EBP[[1]]$y 
      
      list.results$LME_log$est.mean[i,] <- EBP[[1]][[1]]$mean
      list.results$LME_log$est.gini[i,] <- EBP[[1]][[1]]$gini
      list.results$LME_log$est.hcr [i,] <- EBP[[1]][[1]]$hcr 
      list.results$LME_log$est.quant10[i,] <- EBP[[1]][[1]]$quant10
      list.results$LME_log$est.quant25[i,] <- EBP[[1]][[1]]$quant25
      list.results$LME_log$est.quant50[i,] <- EBP[[1]][[1]]$quant50
      list.results$LME_log$est.quant75[i,] <- EBP[[1]][[1]]$quant75
      list.results$LME_log$est.quant90[i,] <- EBP[[1]][[1]]$quant90
      list.results$LME_log$est.pgap[i,] <- EBP[[1]][[1]]$pgap
      list.results$LME_log$est.qsr[i,] <- EBP[[1]][[1]]$qsr
      
      if (MSE==TRUE){
        list.results$LME_log$MSE.mean[i,] <- EBP$MSE[,"mean"]
        list.results$LME_log$MSE.gini[i,] <- EBP$MSE[,"gini"]
        list.results$LME_log$MSE.hcr[i,] <- EBP$MSE[,"hcr"]
        list.results$LME_log$MSE.quant10[i,] <- EBP$MSE[,"quant10"]
        list.results$LME_log$MSE.quant25[i,] <- EBP$MSE[,"quant25"]
        list.results$LME_log$MSE.quant50[i,] <- EBP$MSE[,"quant50"]
        list.results$LME_log$MSE.quant75[i,] <- EBP$MSE[,"quant75"]
        list.results$LME_log$MSE.quant90[i,] <- EBP$MSE[,"quant90"]
        list.results$LME_log$MSE.pgap[i,] <- EBP$MSE[,"pgap"]
        list.results$LME_log$MSE.qsr[i,] <- EBP$MSE[,"qsr"]
      }
    }
    
    #==============================================================================
    # Apply lm boxcox
    
    if ("LME_box" %in% method) {
      source("MSE_function_final_rueck.R")
      EBP <- mse(formulaFix=formel,  formula=formula, population=pop_data, sample=samp_data, L=l, threshold=threshold, 
                 MSE=MSE, B=b, method="LME", transformation="box")
      
      LME_EBP_pop_box <-EBP[[1]]$y 
      
      list.results$LME_box$est.mean[i,] <- EBP[[1]][[1]]$mean
      list.results$LME_box$est.gini[i,] <- EBP[[1]][[1]]$gini
      list.results$LME_box$est.hcr [i,] <- EBP[[1]][[1]]$hcr 
      list.results$LME_box$est.quant10[i,] <- EBP[[1]][[1]]$quant10
      list.results$LME_box$est.quant25[i,] <- EBP[[1]][[1]]$quant25
      list.results$LME_box$est.quant50[i,] <- EBP[[1]][[1]]$quant50
      list.results$LME_box$est.quant75[i,] <- EBP[[1]][[1]]$quant75
      list.results$LME_box$est.quant90[i,] <- EBP[[1]][[1]]$quant90
      list.results$LME_box$est.pgap[i,] <- EBP[[1]][[1]]$pgap
      list.results$LME_box$est.qsr[i,] <- EBP[[1]][[1]]$qsr
      if (MSE==FALSE){LME_Lambda_Sim[i] <- EBP$lambda}
      
      if (MSE==TRUE){
        list.results$LME_box$MSE.mean[i,] <- EBP$MSE[,"mean"]
        list.results$LME_box$MSE.gini[i,] <- EBP$MSE[,"gini"]
        list.results$LME_box$MSE.hcr[i,] <- EBP$MSE[,"hcr"]
        list.results$LME_box$MSE.quant10[i,] <- EBP$MSE[,"quant10"]
        list.results$LME_box$MSE.quant25[i,] <- EBP$MSE[,"quant25"]
        list.results$LME_box$MSE.quant50[i,] <- EBP$MSE[,"quant50"]
        list.results$LME_box$MSE.quant75[i,] <- EBP$MSE[,"quant75"]
        list.results$LME_box$MSE.quant90[i,] <- EBP$MSE[,"quant90"]
        list.results$LME_box$MSE.pgap[i,] <- EBP$MSE[,"pgap"]
        list.results$LME_box$MSE.qsr[i,] <- EBP$MSE[,"qsr"]
        
        LME_Lambda_Sim[i] <- EBP[[3]]$lambda
      }
    }
    
    #  print(paste("=========SimIt=======:", i, "of", simruns))
  }
  
  res <-list(results = list.results , LME_Lambda = LME_Lambda_Sim, EM_Lambda = EM_Lambda_Sim, SEM_Lambda = SEM_Lambda_Sim,
             LME_Prediction = LME_EBP_pop, LME_Prediction_log = LME_EBP_pop_log, LME_Prediction_box = LME_EBP_pop_box,
             EM_Prediction = EM_EBP_pop, EM_Prediction_log = EM_EBP_pop_log, EM_Prediction_box = EM_EBP_pop_box,
             SEM_Prediction = SEM_EBP_pop, SEM_Prediction_log = SEM_EBP_pop_log, SEM_Prediction_box = SEM_EBP_pop_box,
             MID_Prediction = MID_EBP_pop, 
             trueY = true_y, itEM=itEM, itEMLog=itEMLog, itEMBox=itEMBox)
  
  return(res)
  
}


