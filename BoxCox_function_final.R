#Erstmal nur die tranformation testen

# erstmal daten laden
#load("H:/Doktor2/R-Code/Data_transformation.RData")  
#setwd("H:/Doktor2/R-Code/Functions/Paul")

# Formel spezifizieren
# Nur fixe effects
#formel <- as.formula(y~x)
#l <- 5
#formula <- "y ~ x + (1|clusterid)"

boxcox <- function(dat, lambda=NULL, m=NULL, inverse, formula=as.formula("y~x")){

lambda_new <- NULL

#calc_keepPop=function(dat)
#{
#  attr(dat, "pop")=dat
#  dat
#}


#trunc = function(DAT, th)
#{
#  fil = DAT$y < th 
#  DAT$y[fil] = th  
#  return(DAT)
#}

#sampler=function(DAT){
#  smp=as.data.frame(matrix(nrow=sum(sample_size) , ncol=ncol(DAT)))
#  brd=append(0,cumsum(sample_size))
#  for(i in 1:Domains){
#    smp[((brd[i]+1):brd[i+1]),]=(DAT[DAT$idD==i,])[sample(1:sum(DAT$idD==i),size=sample_size[i]),]
#  }
#  attr(smp,"pop")=DAT
#  colnames(smp)=colnames(DAT)
#  return(smp)
#}



geometric.mean<-function (x) #for RMLE in the parameter estimation
{
  exp(mean(log(x)))
}

#

#Transformation Box-Cox 

Box = function(l, y, inv=FALSE, m=NULL) #Box-Cox transformation (lambda=l)
{
  if(!inv)
  {
    if(is.null(m))
    {
      m = 0
    }
    if((s=min(y))<=0) #with shift(=m) parameter for making data positive (>0)
    {
      s = abs(s)+1
    }
    else
    {
      s=0
    }
    m=m+s
    
    if(abs(l)<=1e-12) #case lambda=0
    {
      y = log(y+m)
    }
    else
    {
      y = ((y+m)^l-1)/l
    }
  }
  else
  {
    if(is.null(m)) #inverse transformation
    {
      m = 0
    }
    if(abs(l)<=1e-12) #case lambda=0
    {
      y = exp(y) - m
    }
    else
    {
      y = (l*y+1)^(1/l)-m
    }
  }
  return(list(y = y, m = m)) #return of transformed data and shift (overwriten y)
}


sd_box = function(y, l, m)
{
  if((m=min(y))<=0)
  {
    y=y-m
    y=y+1
  }
  
  gm=geometric.mean(y)
  if(abs(l)>1e-12)
  {
    y=(y^l-1)/(l*((gm)^(l-1)))
  }
  else
  {
    y=gm*log(y)
  }
  return(y)
}








# Estimation cvm

#to_opt_Box <- function(l, y, dat,m){
#  dat$y <- Box(l,y,F,m)$y
#  mod <- lme(y~x, random=~1|idD,method="REML",data = dat)
  #beta schätung unter l fix
#  res=residuals(mod,level=0,type="pearson")
#  step.length=10^-4
#  eval.probs=seq(0,1, by=step.length)
#  eval.points = qnorm(eval.probs,mean=mean(res), sd=sd(res))
#  test.probs = ecdf(res)(eval.points)
#  difs=eval.probs - test.probs # abstand cdf norm zu ecd
  #print(select_john(z = 1, res))
#  result = sum((difs)^2)#cramer von mises
#  return(result)
#}



# Estimation REML

opt_reml_box <- function(l,y,dat,m)
{
  dat$y = sd_box(y=y, l = l, m=m)
  mod <- lme(formula, random=~1|clusterid,method="REML",data = dat, control = lmeControl(opt = "optim"))
  if(logLik(mod)>0)
  {
    return(99999)
  }
  return(-logLik(mod))
}

if (is.data.frame(dat)==TRUE) {
  if (inverse==FALSE){
    l_box <- optimise(f = opt_reml_box, interval = c(-1,2), 
                    y= dat$y,dat = dat, m = NULL)$minimum
    trans_dat <- Box(l_box, y = dat$y, inv=FALSE)$y
    m <- Box(l_box, y = dat$y, inv=FALSE)$m
    lambda_new <- l_box
  }

  if (inverse==TRUE){
    trans_dat <- Box(lambda, y=dat$y , inv = TRUE, m=m)$y
  }
}

if (is.data.frame(dat)==FALSE) {
  if (inverse==FALSE){
    trans_dat <- Box(lambda, y = dat, inv=FALSE)$y
    m <- Box(lambda, y = dat, inv=FALSE)$m
    
  }
  
  if (inverse==TRUE){
    trans_dat <- Box(lambda, y=dat , inv = TRUE, m=m)$y
  }
}

res <- list(trans_dat, m, lambda_new)
return(res)

}





