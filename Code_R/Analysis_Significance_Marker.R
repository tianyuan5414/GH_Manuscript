#' This script aims at adding markers (*) according to the significance level and criteria
sigMarkter <- function(tempSig, gradient) {
  
  #Significant level 1, marked by *
  if ((as.numeric(tempSig) < gradient[1]) & 
      (as.numeric(tempSig) > gradient[2])) {
    
    return('*')
    
    #Significance level 2, marked by **  
  } else if ((as.numeric(tempSig) < gradient[2]) & 
             (as.numeric(tempSig) > gradient[3])) {
    
    return('**')
    
    #Significance level 3, marked by ***
  } else if (as.numeric(tempSig) < gradient[3]) {
    
    return('***')
    
  } else {
    
    #Not significant, marked by nothing
    return('')
    
  }
}
