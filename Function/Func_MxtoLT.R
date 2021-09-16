Fun_MxtoLT <- function(x, Mx, sex = "M", ax = NULL, outcome){
  m <- length(x)
  n <- c(diff(x), 1)
  
  if(is.null(ax)){
    ax <- rep(0, m)
    if(x[1] != 0 | x[2] != 1){
      ax <- n / 2
      ax[m] <- 1 / Mx[m]
    }else{    
      if(sex == "F"){
        if(Mx[1] >= 0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800 * Mx[1]
        }
      }
      if(sex == "M"){
        if(Mx[1] >= 0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684 * Mx[1]
        }
      }
      ax[-1] <- n[-1] / 2
      ax[m] <- 1 / Mx[m]
    }
  }
  qx  <- n * Mx / (1 + (n - ax) * Mx)
  qx[m] <- 1
  px  <- 1 - qx
  lx  <- cumprod(c(1, px))
  dx  <- -diff(lx)
  Lx  <- n * lx[-1] + ax * dx
  lx <- lx[-(m + 1)]
  Lx[m] <- lx[m] / Mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx / lx
  
  if(outcome == "lx"){
    return(lx)
    
  } else if(outcome == "ex"){
    return(ex)
    
  } else {
    return(qx)
  }
  
}