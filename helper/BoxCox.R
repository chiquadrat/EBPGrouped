#Author: Paul Walter
#Mail: paul.walter@fu-berlin.de

boxcox <- function(dat, lambda=NULL, m=NULL, inverse, formula=as.formula("y~x")){

lambda_new <- NULL



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





