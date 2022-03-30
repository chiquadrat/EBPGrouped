#Author: Paul Walter
#Mail: paul.walter@fu-berlin.de

SEM_lme <- function(data, classes, formula, burnin = 2, samples = 5) {


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
    means[negInf] <- sapply(intervals[negInf], function(x) (x[2] +  (x[2]-
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
  std.error_test <- unname((summary(regclass)$coefficients)[,2])


  resulty <- matrix(ncol = c(burnin + samples), nrow = length(yclass))
  resultcoef <- matrix(ncol = c(burnin + samples), nrow = length(regclass@beta))
  result_ranef <- matrix(ncol = c(burnin + samples), nrow = length(unname(unlist(ranef(regclass)))))
  result_sigmae<-vector(mode = "numeric", length = burnin+samples)
  result_sigmau<-vector(mode = "numeric", length = burnin+samples)
  result_std.error <- matrix(ncol = c(burnin+samples), nrow = length(unname((summary(regclass)$coefficients)[,2])))

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
    resultcoef[,j] <- regclass@beta
    result_ranef[,j]<-unname(unlist(ranef(regclass)))
    result_sigmae[j]<- sigmahat
    result_sigmau[j]<-as.data.frame(VarCorr(regclass))$sdcor[1]
    resulty[,j] <- data$pseudoy
    result_std.error[,j] <- unname((summary(regclass)$coefficients)[,2])
    print(paste("Iteration:", j, "of", burnin + samples))
  }

  est <- list(data=data, resultY = resulty, coef = resultcoef, ranef = result_ranef, sigmae = result_sigmae,
              sigmau = result_sigmau, midreg = regclass0@beta, std.error = result_std.error, std.errorTest=std.error_test )
  print("SEM done")
  return(est)
}




