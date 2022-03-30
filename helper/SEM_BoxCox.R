#Author: Paul Walter
#Mail: paul.walter@fu-berlin.de

SEM_lme_BoxCox <- function(data, classes, formula, burnin = 2, samples = 5) {

  data2 <- data
  
  yclass <- data$yclass
  yclassl <- data$yclass
  levels(yclassl) <- 1:length(levels(yclassl))
  yclassl <- as.numeric(as.vector(yclassl))
  data$yclassl <- yclassl
  
  
  formula <- as.formula(gsub(".*~","pseudoy~",formula))
  
  
  intervals <- vector("list", length(classes) -  1)
  
  
  for (i in seq(length = length(classes) - 1)) { 
    intervals[[i]] <- c(classes[i], classes[i +1])
  }
  
  means <- sapply(intervals, mean) 
  widths <- sapply(intervals, function(x) x[2] - x[1]) 
  meanWidth <- mean(widths[!is.infinite(widths)]) 
  negInf <- is.infinite(means) & means < 0 
  
  if (any(negInf)) {
    means[negInf] <- sapply(intervals[negInf], function(x) (x[2] -  (x[2]-
                                                                       meanWidth))/2)
  }
  posInf <- is.infinite(means) & means > 0 
  if (any(posInf)) {
    means[posInf] <- sapply(intervals[posInf], function(x) (x[1]+(x[1] + 
                                                                    meanWidth))/2)
  }
  
  
  yclassmeans<-means
  
  levels(yclass) <- yclassmeans
  data$pseudoy <- as.numeric(as.vector(yclass))
  
  formelBoxCox <- as.character(formula)
  formelBoxCox <- sub("\\s+\\+\\s+\\(.*", "", formelBoxCox)
  formelBoxCox <- gsub("pseudoy", "y", formelBoxCox)
  formelBoxCox <- as.formula(formelBoxCox)
  
  data$yold <- data$y
  data$y <- data$pseudoy
  source("helper/BoxCox.R")
  BoxCox <- boxcox(dat=data, inverse=FALSE, formula = formelBoxCox)
  BoxCox_y <- BoxCox[[1]]
  BoxCox_m <- BoxCox[[2]]
  BoxCox_lambda <- BoxCox[[3]]
  data$pseudoy <- BoxCox_y
  

  source("helper/BoxCox.R")
  BoxCoxClasses <- boxcox(dat=classes,lambda = BoxCox_lambda, inverse=FALSE)
  classesBox <- BoxCoxClasses[[1]]
  classesM <- BoxCoxClasses[[2]]
  
  regclass0 <- lmer(formula,data=data )
  regclass <- lmer(formula,data=data )
  
  it_lambda <- (burnin+samples)*5
  
  resulty <- matrix(ncol = c(it_lambda), nrow = length(yclass))
  resultcoef <- matrix(ncol = c(it_lambda), nrow = length(regclass@beta))
  result_ranef <- matrix(ncol = c(it_lambda), nrow = length(unname(unlist(ranef(regclass)))))
  result_sigmae<-vector(mode = "numeric", length = it_lambda)
  result_sigmau<-vector(mode = "numeric", length = it_lambda)
  result_lambda <- vector(mode = "numeric", length = it_lambda)
  result_m <- vector(mode = "numeric", length = it_lambda)
  
  for (j in 1:(it_lambda)) {
    data$predict <- predict(regclass,data)
    sigmahat <- sigma(regclass)
    for (i in 1:(length(classesBox) - 1)) {
      if (nrow(data[data$yclassl==i,])!=0) { 
      mean <- data$predict[data$yclassl==i]
      pseudoy <- rtruncnorm(length(mean), a=classesBox[i], b=classesBox[i+1], mean=mean, sd=sigmahat )
      data$pseudoy[data$yclassl==i] <- pseudoy
      }
    }
    
    result_lambda[j] <- BoxCox_lambda
    result_m[j] <- BoxCox_m
    
    data$y <- data$pseudoy
    rueck <- boxcox(dat=data, m=BoxCox_m, lambda = BoxCox_lambda, inverse = T)
    data$y <- rueck[[1]]
    resulty[,j] <- rueck[[1]]
    
    source("helper/BoxCox.R")
    BoxCox <- boxcox(dat=data, inverse=FALSE, formula = formelBoxCox)
    BoxCox_y <- BoxCox[[1]]
    BoxCox_m <- BoxCox[[2]]
    BoxCox_lambda <- BoxCox[[3]]
    data$pseudoy <- BoxCox_y
    
    source("helper/BoxCox.R")
    BoxCoxClasses <- boxcox(dat=classes,lambda = BoxCox_lambda,  inverse=FALSE)
    classesBox <- BoxCoxClasses[[1]]
    classesM <- BoxCoxClasses[[2]]
    
    regclass=lmer(formula,data=data )
    
    resultcoef[,j] <- regclass@beta
    result_ranef[,j]<-unname(unlist(ranef(regclass)))
    result_sigmae[j]<- sigmahat
    result_sigmau[j]<-as.data.frame(VarCorr(regclass))$sdcor[1]                  
  }
  
  lambda <- mean(result_lambda[-c(1:(burnin*5))])
  
  source("helper/BoxCox.R")
  BoxCoxClasses <- boxcox(dat=classes,lambda = lambda,  inverse=FALSE)
  classesBox <- BoxCoxClasses[[1]]
  classesM <- BoxCoxClasses[[2]]
  
  classes <- classesBox
  data <- data2
  
  yclass <- data$yclass
  yclassl <- data$yclass
  levels(yclassl) <- 1:length(levels(yclassl))
  yclassl <- as.numeric(as.vector(yclassl))
  data$yclassl <- yclassl
  
  intervals <- vector("list", length(classes) -  1)
  
  for (i in seq(length = length(classes) - 1)) { 
    intervals[[i]] <- c(classes[i], classes[i +1])
  }
  
  means <- sapply(intervals, mean) 
  widths <- sapply(intervals, function(x) x[2] - x[1]) 
  meanWidth <- mean(widths[!is.infinite(widths)]) 
  negInf <- is.infinite(means) & means < 0 
  
  if (any(negInf)) {
    means[negInf] <- sapply(intervals[negInf], function(x) (x[2] -  (x[2]-
                                                                       meanWidth))/2)
  }
  posInf <- is.infinite(means) & means > 0 
  if (any(posInf)) {
    means[posInf] <- sapply(intervals[posInf], function(x) (x[1]+(x[1] + 
                                                                    meanWidth))/2)
  }
  
  
  yclassmeans<-means
  
  levels(yclass) <- yclassmeans
  data$pseudoy <- as.numeric(as.vector(yclass))
  
  regclass0 <- lmer(formula,data=data )
  regclass <- lmer(formula,data=data )
  
  resulty_step2 <- matrix(ncol = c(burnin + samples), nrow = length(yclass))
  resultcoef_step2 <- matrix(ncol = c(burnin + samples), nrow = length(regclass@beta))
  result_ranef_step2 <- matrix(ncol = c(burnin + samples), nrow = length(unname(unlist(ranef(regclass)))))
  result_sigmae_step2 <- vector(mode = "numeric", length = burnin+samples)
  result_sigmau_step2 <- vector(mode = "numeric", length = burnin+samples)
  result_std.error_step2 <- matrix(ncol = c(burnin+samples), nrow = length(unname((summary(regclass)$coefficients)[,2])))
  
  for (j in 1:(burnin + samples)) {
    data$predict <- predict(regclass,data)
    sigmahat <- sigma(regclass)
    for (i in 1:(length(classes) - 1)) {
      if (nrow(data[data$yclassl==i,])!=0) { 
        mean <- data$predict[data$yclassl==i]
        pseudoy <- rtruncnorm(length(mean), a=classes[i], b=classes[i+1], mean=mean, sd=sigmahat )
        data$pseudoy[data$yclassl==i] <- pseudoy
      }
    }
    regclass=lmer(formula,data=data )
    resultcoef_step2[,j] <- regclass@beta
    result_ranef_step2[,j]<-unname(unlist(ranef(regclass)))
    result_sigmae_step2[j]<- sigmahat
    result_sigmau_step2[j]<-as.data.frame(VarCorr(regclass))$sdcor[1]                  
    resulty[,j] <- data$pseudoy
    result_std.error_step2[,j] <- unname((summary(regclass)$coefficients)[,2])
  }
  
  est <- list(data=data, resultY = resulty, coef = resultcoef_step2, ranef = result_ranef_step2, 
              sigmae = result_sigmae_step2, sigmau = result_sigmau_step2, midreg = regclass0@beta, 
              lambda = result_lambda, m=result_m, coef_pre = resultcoef, ranef_pre = result_ranef, 
              sigmae_pre = result_sigmae, sigmau_pre = result_sigmau, std.error = result_std.error_step2)
  print("SEM Box Cox done")
  return(est)
}




