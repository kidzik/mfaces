igcv.wrapper <-
function(Y, X, P, N2, lower, upper, search.length, 
                         search.length0 = search.length){
  
  Rho <- exp(seq(lower,upper,length=search.length))
  List <- seq(0.01,0.99,length=search.length0) 
  fit_list <- list(length=length(List))
  # Time <- proc.time()
  BtB <- crossprod(X)
  EigenB <- eigen(BtB)
  EigenB$values <- sapply(EigenB$values, function(x) max(x,0.00000001))
  B_inv_half <- EigenB$vectors%*%diag(sqrt(1/EigenB$values))%*%t(EigenB$vectors)
  
  P1 <- B_inv_half%*%P[[1]]%*%B_inv_half
  P2 <- B_inv_half%*%P[[2]]%*%B_inv_half
  
  B <- as.matrix(X %*% B_inv_half)
  f = crossprod(B,Y) 
  c2 = ncol(BtB)
  
  mat_list <- list()
  for(i in 1:length(N2)){
    
    seq = (sum(N2[1:i])-N2[i]+1):(sum(N2[1:i]))
    Bi = matrix(B[seq,], nrow = length(seq))  
    fi = crossprod(Bi,Y[seq])
    Li = crossprod(Bi)
    
    LList <- list()
    LList[[1]] = fi
    LList[[2]] = Li
    mat_list[[i]] = LList
    
  }
  # print(proc.time()-Time)
  # Time <- proc.time()
  for(i in 1:length(List)){
    w = List[i]
    fit_list[[i]] <- igcv.criteria(B,Y,mat_list,f,P1,P2,N2,Rho,w,c2)
  }
  # print(proc.time()-Time)
  index <- which.min(sapply(fit_list,function(x) x$igcv))
  res <- fit_list[[index]]
  res$lambda <- c(List[index]*res$rho, (1-List[index])*res$rho)

  return(res)
}
