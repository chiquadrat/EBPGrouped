mse = function(formulaFix, formula, population , sample, L=5, threshold=0.5, MSE=FALSE , B=5,
               burnin=5, samples=20, classes, method="LME", terminate = 0.01, transformation ="none"){
  
  
  # General Variable Definitions ---------------------------------
   
  uni_domains <- population$clusterid
  rearrange_universe = order(uni_domains)
  
  smp_domains <- sample$clusterid
  rearrange_sample = order(smp_domains)
  
  universe = population[rearrange_universe,]
  uni_domains = uni_domains[rearrange_universe]
  
  
  sample = sample[rearrange_sample,]
  smp_domains = smp_domains[rearrange_sample]
  
  
  N_uni = length(uni_domains)
  N_smp = length(smp_domains)
  N_unobs= N_uni-N_smp
  N_dom_smp = length(unique(smp_domains))
  N_dom_uni = length(unique(uni_domains))
  N_dom_unobs = N_dom_uni - N_dom_smp
  n_smp = as.vector(table(smp_domains))
  n_uni = as.vector(table(uni_domains))
  
  obs_dom = uni_domains %in% unique(smp_domains)
  dist_obs_dom = sort(unique(uni_domains)) %in% unique(smp_domains)
  
  # test
  tmp = formulaFix
  tmp[2] = formula(uni_domains~1)[2]
  X_uni = model.matrix(tmp, data.frame(uni_domains, universe))
  rm(tmp)
  X_smp = model.matrix(formulaFix,sample)
  Y_smp = as.matrix(sample[paste(formulaFix[2])])
  
  
  # PovertyFunctions ---------------------------------------------------------
  
  hcr_function = function(y){
    mean(y<threshold)
  }
  pgap_function = function(y){
    mean((y<threshold)*(threshold-y)/threshold)
  }
  qsr_function = function(y){
    sum(y[(y>quantile(y,0.8))])/sum(y[(y<quantile(y,0.2))])
  }
  fast_mean = function(X)
  {
    .Internal(mean(X))
  }
  gini_function = function(X)
  {
    n <- length(X)
    X <- sort(X)
    G <- sum(X*1L:n)
    G <- 2*G/sum(X) - (n + 1L)
    return(G/n)
  }
  
  # Point Estimation Function ------------------------------------
  point_estim <- function(Y_smp, classify) {
  
  # Estimation of the mixed linear model
  sampleEBP <- sample
  sampleEBP$y <- unname(Y_smp) 
  
  if (classify==TRUE & method!="LME") {
    sampleEBP$yclass <- cut(Y_smp, classes)
  }
  
  
  if (transformation=="none") {
    
    if (method=="MID") {    
      
      source("MID_lme_function.R")
      MID <- MID_lme(data = sampleEBP, classes, formula)
      coef <- MID$coef
      ranef <- MID$ranef
      sigmae <- MID$sigmae
      sigmau <- MID$sigmau
      m <- NULL
      lambda <- NULL
      itEM <- NULL
    }
    
    
    if (method=="LME") {    
      lmreg <- lmer(formula, data=sampleEBP)  
      coef <- lmreg@beta
      ranef <- unname(unlist(ranef(lmreg)))
      sigmae <- sigma(lmreg)
      sigmau <- unname(as.data.frame(VarCorr(lmreg))$sdcor[1])
      m <- NULL
      lambda <- NULL
      itEM <- NULL
    }
  
    if (method=="SEM") {
      
      source("SEM_lme_function.R")
      SEM <- SEM_lme(data=sampleEBP, classes=classes, formula, burnin, samples)
      coef <- apply(SEM$coef[,-c(1:burnin)],1,mean)
      ranef <- apply(SEM$ranef[,-c(1:burnin)],1,mean)
      sigmae <- mean(SEM$sigmae[-c(1:burnin)])
      sigmau <- mean(SEM$sigmau[-c(1:burnin)])
      m <- NULL
      lambda <- NULL
      itEM <- NULL
    }
  
    if (method=="EM") {
      source("EM_lme_function.R")
      EM <- EM_lme(data=sampleEBP, classes=classes, formula, terminate = terminate)
      coef <- EM$coef[,ncol(EM$coef)]
      ranef <- EM$ranef[,ncol(EM$ranef)]
      sigmae <- EM$sigmae[,ncol(EM$sigmae)]
      sigmau <- EM$sigmau[,ncol(EM$sigmau)]
      m <- NULL
      lambda <- NULL
      itEM <- EM$iteration
    }
  }
  
  if (transformation=="log") {
    if (method=="LME") {    
      
      # Transform data
      source("Log_function_final.R")
      daten <- log_function(y = sampleEBP$y)
      sampleEBP$y <- daten[[1]]
      m <- daten[[2]]
      lambda <- NULL
      itEM <- NULL
      
      lmreg <- lmer(formula, data=sampleEBP)  
      coef <- lmreg@beta
      ranef <- unname(unlist(ranef(lmreg)))
      sigmae <- sigma(lmreg)
      sigmau <- unname(as.data.frame(VarCorr(lmreg))$sdcor[1]) }
    
    if (method=="SEM") {
      # Transform classes
      source("Log_function_final.R")
      daten <- log_function(y = classes)
      classesLog <- daten[[1]]
      m <- daten[[2]]
      lambda <- NULL
      itEM <- NULL
            
      source("SEM_lme_function.R")
      SEM <- SEM_lme(data=sampleEBP, classesLog, formula, burnin, samples)
      coef <- apply(SEM$coef[,-c(1:burnin)],1,mean)
      ranef <- apply(SEM$ranef[,-c(1:burnin)],1,mean)
      sigmae <- mean(SEM$sigmae[-c(1:burnin)])
      sigmau <- mean(SEM$sigmau[-c(1:burnin)])
    }
    
    if (method=="EM") {
      # Transform classes
      source("Log_function_final.R")
      daten <- log_function(y = classes)
      classesLog <- daten[[1]]
      m <- daten[[2]]
      lambda <- NULL
            
      source("EM_lme_function.R")
      EM <- EM_lme(data=sampleEBP, classesLog, formula, terminate = terminate)
      coef <- EM$coef[,ncol(EM$coef)]
      ranef <- EM$ranef[,ncol(EM$ranef)]
      sigmae <- EM$sigmae[,ncol(EM$sigmae)]
      sigmau <- EM$sigmau[,ncol(EM$sigmau)]
      itEM <- EM$iteration
    }
  }
  
  if (transformation=="box") {
    if (method=="LME") {    
      source("BoxCox_function_final.R")
      sampBox <- boxcox(dat=sampleEBP, inverse = FALSE, formula = formulaFix)
      #samp_dataBox <- samp_data
      sampleEBP$y <- sampBox[[1]]
      m <- sampBox[[2]]
      lambda <- sampBox[[3]]
      itEM <- NULL
      
      lmreg <- lmer(formula, data=sampleEBP)
      
      coef <- lmreg@beta
      ranef <- unname(unlist(ranef(lmreg)))
      sigmae <- sigma(lmreg)
      sigmau <- unname(as.data.frame(VarCorr(lmreg))$sdcor[1])
      
    }
    
    if (method=="SEM") {
      source("SEM_lme_function_BoxCox.R")
      SEM <- SEM_lme_BoxCox(data=sampleEBP, classes, formula, burnin, samples)
      coef <- apply(SEM$coef[,-c(1:burnin)],1,mean)
      ranef <- apply(SEM$ranef[,-c(1:burnin)],1,mean)
      sigmae <- mean(SEM$sigmae[-c(1:burnin)])
      sigmau <- mean(SEM$sigmau[-c(1:burnin)])
      m <- SEM$m[1]
      lambda <- mean(SEM$lambda[-c(1:(burnin*5))])
      itEM <- NULL
    }
    
    if (method=="EM") {
      source("EM_lme_function_BoxCox.R")
      EM <- EM_lme_BoxCox(data=sampleEBP, classes, formula, terminate = terminate)
      coef <- EM$coef[,ncol(EM$coef)]
      ranef <- EM$ranef[,ncol(EM$ranef)]
      sigmae <- EM$sigmae[,ncol(EM$sigmae)]
      sigmau <- EM$sigmau[,ncol(EM$sigmau)]
      m <- EM$m[1]
      lambda <- EM$lambda[length(EM$lambda)]
      itEM <- EM$iteration
      
      
    }
  }
  
  
  betas=coef
  rand_eff=rep(0, length(unique(uni_domains)))            
  rand_eff[dist_obs_dom]=ranef
  sigmae2est=sigmae^2    
  sigmau2est=sigmau^2
  
  gamma = sigmau2est / (sigmau2est + sigmae2est / n_smp) 
  sigmav2est = sigmau2est * (1-gamma)
  
  rand_eff_uni = rep(rand_eff, n_uni)  
  #rand_eff_smp = rep(rand_eff, n_smp) #NEU
  
  mu = X_uni %*% betas + rand_eff_uni
  
  quant10s = quant25s = quant50s = quant75s = quant90s = ginis = qsrs = pgaps = hcrs = means = matrix(nrow=N_dom_uni, ncol=L) 
  
  # Save EBP prediction
  y_EBP <- vector("list",1)
  
  # for-loop calculating the pseudo populations
  for(l in 1 : L)
  {
    eps = rnorm(N_uni, 0, sqrt(sigmae2est))
    
    vu=vector(length=N_uni)
    #browser()
    vu[!obs_dom] = rep(rnorm(N_dom_unobs, 0, sqrt(sigmau2est)),n_uni[!dist_obs_dom])
    vu[obs_dom]  = rep(rnorm(rep(1, N_dom_smp), 0, sqrt(sigmav2est)), n_uni[dist_obs_dom])
    y_pred = mu + eps + vu
    
    if (transformation=="log") {
      source("Log_function_final.R")
      y_pred <- log_function(y=y_pred, inv=T, m=m)[[1]]
    }
    
    if (transformation=="box") {
      y_pre_vec <- y_pred[,1]
      source("BoxCox_function_final.R")
      y_pre_vec <- boxcox(dat=y_pre_vec, m=m, lambda = lambda,  inv=T)[[1]]
      y_pred[,1]<-y_pre_vec
    }
    
    y_pred[!is.finite(y_pred)] = 1.1
    
    y_EBP[[1]] <- y_pred
    
    quant10s[,l] = tapply(y_pred, uni_domains, quantile, probs=0.1)
    quant25s[,l] = tapply(y_pred, uni_domains, quantile, probs=0.25)
    quant50s[,l] = tapply(y_pred, uni_domains, quantile, probs=0.50)
    quant75s[,l] = tapply(y_pred, uni_domains, quantile, probs=0.75)
    quant90s[,l] = tapply(y_pred, uni_domains, quantile, probs=0.9)
    ginis[,l] = tapply(y_pred, uni_domains, gini_function)
    means[,l] = tapply(y_pred, uni_domains, fast_mean)
    hcrs[,l] = tapply(y_pred, uni_domains, hcr_function)
    qsrs[,l] =  tapply(y_pred, uni_domains, qsr_function)
    pgaps[,l] =  tapply(y_pred, uni_domains, pgap_function)
  }
  
  #plot(density(quant90s))
  
  # point estimations
  results = list(data.frame(
    dom = unique(uni_domains),
    quant10 = rowMeans(quant10s),
    quant25 = rowMeans(quant25s),
    quant50 = rowMeans(quant50s),
    quant75 = rowMeans(quant75s),
    quant90 = rowMeans(quant90s),
    gini = rowMeans(ginis),
    mean = rowMeans(means),
    hcr = rowMeans(hcrs),
    qsr = rowMeans(qsrs),
    pgap = rowMeans(pgaps)
  ),  y = y_EBP)
  
  return(list(POV=results,
              lambda=lambda,
              m=m,
              mu=mu,
              betas=betas,
              rand_eff=rand_eff,
              rand_eff_uni=rand_eff_uni,
              #rand_eff_smp=rand_eff_smp,
              sigmae2est=sigmae2est,
              sigmau2est=sigmau2est,
              sigmav2est=sigmav2est,
              it = itEM )
  )
  
  }
  
 # bla <- point_estim()

# MSE Estim Function--------------------------------------------
  mse_estim = function(mu_uni_, mu_smp_,sigmae2est_,sigmau2est_,sigmav2est_ ,m_ , lambda_) 
  {
    eps   =  vector(length = N_uni )
    eps[obs_dom]  =  rnorm(sum(obs_dom), 0, sqrt(sigmae2est_))
    eps[!obs_dom] =  rnorm(sum(!obs_dom), 0, sqrt(sigmae2est_+sigmau2est_))
    vu_tmp=rnorm(N_dom_uni, 0, sqrt(sigmau2est_))
    vu_uni  = rep(vu_tmp, n_uni)
    vu_smp  = rep(vu_tmp[dist_obs_dom], n_smp)
    rm(vu_tmp)
    Y_uni_b = mu_uni_ + eps + vu_uni
    eps = rnorm(N_smp, 0, sqrt(sigmae2est_))

    Y_smp_b = mu_smp_ + eps + vu_smp
    #Y_smp_b = transformation(y=Y_smp_b, l=optpar_, inv=T, m=par_m_)$y
    
    if (transformation=="log") {
      source("Log_function_final.R")
      Y_smp_b <- log_function(y=Y_smp_b, inv=T, m=m_)[[1]]
    }
    
    if (transformation=="box") {
      y_pre_vec <- Y_smp_b[,1]
      source("BoxCox_function_final.R")
      y_pre_vec <- boxcox(dat=y_pre_vec, m=m_, lambda = lambda_,  inv=T)[[1]]
      Y_smp_b[,1]<-y_pre_vec
    }
    
    Y_smp_b[!is.finite(Y_smp_b)]=1.1
    Y_smp_b[Y_smp_b<=1] = 1.1

    estims = as.matrix(point_estim(Y_smp = Y_smp_b , classify = "TRUE")[[1]][[1]][,-1])
    
    #Y_uni_b = transformation(y=Y_uni_b, l=optpar_, inv=T, m=par_m_)$y
    if (transformation=="log") {
      source("Log_function_final.R")
      Y_uni_b <- log_function(y=Y_uni_b, inv=T, m=m_)[[1]]
    }
    
    if (transformation=="box") {
      y_pre_vec <- Y_uni_b[,1]
      source("BoxCox_function_final.R")
      y_pre_vec <- boxcox(dat=y_pre_vec, m=m_, lambda = lambda_,  inv=T)[[1]]
      Y_uni_b[,1]<-y_pre_vec
    }
    
    Y_uni_b[!is.finite(Y_uni_b)] = 1.1
    Y_uni_b[Y_uni_b<=1] = 1.1
    truth = as.matrix(data.frame(
      quant10 = tapply(Y_uni_b, uni_domains, quantile, probs=0.1),
      quant25 = tapply(Y_uni_b, uni_domains, quantile, probs=0.25),
      quant50 = tapply(Y_uni_b, uni_domains, quantile, probs=0.50),
      quant75 = tapply(Y_uni_b, uni_domains, quantile, probs=0.75),
      quant90 = tapply(Y_uni_b, uni_domains, quantile, probs=0.9),
      gini = tapply(Y_uni_b, uni_domains, gini_function),
      mean = tapply(Y_uni_b, uni_domains, fast_mean),
      hcr =  tapply(Y_uni_b, uni_domains, hcr_function),
      qsr = tapply(Y_uni_b, uni_domains, qsr_function),
      pgap = tapply(Y_uni_b, uni_domains, pgap_function)
    ))

   return((estims-truth)^2)
  }
  

# Calling the functions ----------------------------------------
  
  Pov_ests = point_estim(Y_smp = Y_smp, classify = "FALSE") 

  mu_smp = X_smp %*% Pov_ests$betas #+ Pov_ests$rand_eff_smp # fehlt hier der random effekt?!?!
  mu_uni = X_uni %*% Pov_ests$betas #+ Pov_ests$rand_eff_uni # fehlt hier der random effekt?!?!

  if (MSE)
  {
    mses = replicate(B, mse_estim(
      lambda_         = Pov_ests$lambda ,  
      m_              = Pov_ests$m ,
      mu_uni_         = mu_uni ,
      mu_smp_         = mu_smp ,
      sigmae2est_     = Pov_ests$sigmae2est ,
      sigmau2est_     = Pov_ests$sigmau2est ,
      sigmav2est_     = Pov_ests$sigmav2est
      ))
    mses = apply(mses, c(1,2), fast_mean)
    mses = cbind(unique(uni_domains),mses)
    result = list(POV = Pov_ests[[1]], MSE = mses, Pov_ests[-1])
  } else
  {
    result = Pov_ests
  }

  return(result)

}



