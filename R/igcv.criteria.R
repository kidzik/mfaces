igcv.criteria <-
function(B,Y,mat_list,f,P1,P2,N2,Rho,w,c2){
  S <- w * P1 + (1 - w) * P2
  EigenS <- eigen(S)
  
  s <- EigenS$values
  s[s<0.00000001] <- 0
  U = EigenS$vectors
  tU = t(EigenS$vectors)
  
  ftilde = tU %*% f
  tftilde = t(ftilde)
  
  g <- rep(0, c2)
  G <- matrix(0,c2,c2)
  Li_list <- list()

  # Time <- proc.time()
  for(i in 1:length(N2)){
    
    fi = tU %*% mat_list[[i]][[1]]
    # fi = mat_list[[i]][[1]]
    Li = as.matrix(tU %*% mat_list[[i]][[2]]) %*% U
    # Li = mat_list[[i]][[2]]
    
    g = g + fi*fi
    G = G + Li*(fi%*%tftilde)
    
    Li_list[[i]] = Li
    
  }

  igcv <- function(x){
    d <- 1/(1+x*s)
    ftilde_d <- ftilde*d
    cv0 <- -2*sum(ftilde_d*ftilde)
    cv1 <-  sum(ftilde_d^2)
    cv2 <-  2*sum(d*g)
    cv3 <-  -4*sum(d*(G%*%d))
    cv4 <- sum(unlist(sapply(Li_list,function(x){
      a <- x%*%ftilde_d
      2*sum(a^2*d)
    })))
    cv <- cv0 + cv1 + cv2 + cv3 + cv4
    return(cv)
  }
  
  Length <- length(Rho)
  iGcv <- rep(0,Length)
  # Time <- proc.time()
  for(i in 1:Length) 
    iGcv[i] <- igcv(Rho[i])
  # print(proc.time()-Time)
  index <- which.min(iGcv)
  rho.min <- Rho[index]
  res <- list("rho"=rho.min, "igcv" = iGcv[index])
  return(res)
}
