#==================================================
#
# Function regression on midpoints
#
#==================================================

MID_lme <- function(data, classes, formula) {
  
  samp_data <- data
  intervals <- vector("list", length(classes) -  1)# Liste mit intervalgrenzen
  
  
  for (k in seq(length = length(classes) - 1)) { # Liste mit intervalgrenzen erstellen
    intervals[[k]] <- c(classes[k], classes[k +1])
  }
  
  means <- sapply(intervals, mean) # Klassenmittelwerte
  widths <- sapply(intervals, function(x) x[2] - x[1]) #Klassenbreiten
  meanWidth <- mean(widths[!is.infinite(widths)]) # Durchschnittliche Klassenbreite (ohne unendlich)
  negInf <- is.infinite(means) & means < 0 # Gibt es negative unendliche klassen
  
  if (any(negInf)) {
    means[negInf] <- sapply(intervals[negInf], function(x)(x[2] -  (x[2]-
                                                                      meanWidth))/2)
  }
  posInf <- is.infinite(means) & means > 0 # Gibt es positive unendliche klassen
  if (any(posInf)) {
    means[posInf] <- sapply(intervals[posInf], function(x) (x[1]+(x[1] + # wenn ja, dann Inf ersetzten mit durchschnittlicher klassenbreite
                                                                    meanWidth))/2)
  }
  
  
  yclassmeans<-means
  samp_data$yclassadj <- samp_data$yclass 
  levels(samp_data$yclassadj) <- yclassmeans
  samp_data$ymid <- as.numeric(as.vector(samp_data$yclassadj))
  
  formula_mid <- as.formula(gsub(".*~","ymid~",formula))
  
  lmmidreg <- lmer(formula_mid,data=samp_data)
  
  coef <- lmmidreg@beta
  ranef <- unname(unlist(ranef(lmmidreg)))
  sigmae <- sigma(lmmidreg)
  sigmau <- unname(as.data.frame(VarCorr(lmmidreg))$sdcor[1])
  std.error <- unname((summary(lmmidreg)$coefficients)[,2])
  
  est <- list(coef = coef, ranef = ranef, sigmae = sigmae, sigmau = sigmau, std.error = std.error)
  print("MID done")
  return(est)
  
}