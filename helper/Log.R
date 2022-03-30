#Author: Paul Walter
#Mail: paul.walter@fu-berlin.de

log_function <- function(y, inv = FALSE,  m=0) {
  
  
    if(!inv)
    {
      m = 0
      if(min(y)<1)
      {
        m = 1
        y = y + m
      }
      y=log(y)
    } else
    {
      y=exp(y)-m
    }
    return(list(y=y, m=m))
  
}
  
 


