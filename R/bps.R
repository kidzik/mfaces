bps <-
function(y, x, z, N2, knots=10, p=3, m=2, lower=-3, upper=10, 
                search.length=14,
                knots.option="equally-spaced",lambda=NULL){
  
  
  data <- data.frame("x" = x, "z" = z)
  
  # bs <- rep(bs,2)[1:2]
  knots <- rep(knots,2)[1:2]
  p <- rep(p,2)[1:2]
  m <- rep(m,2)[1:2]
  knots.option <- rep(knots.option,2)[1:2]

  ## x-axis
  x.knots <- construct.knots(data$x,knots[1],knots.option[1],p[1])
  x.List <- pspline.setting(data$x,knots=x.knots,p[1],m[1],
                            type="simple",knots.option=knots.option[1])
  
  # x.object <- s(x,bs=bs[1],k=k[1])
  # x.knots <- NULL
  # if(knots.option[1] == "quantile") x.knots <- select.knots(x,knots = k[1] -3)
  
  # object <- smooth.construct(x.object, 
  # data = data.frame(x = data$x),knots=list(x=x.knots))
  
  Bx <- x.List$B
  Px <- x.List$P
  cx <- nrow(Px)
  
  x.object = list(x.List=x.List,
                  knots = x.knots,
                  knots.option = knots.option[1])
  
  
  ## z-axis
  # z.object <- s(z,bs=bs[2],k=k[2])
  # z.knots <- NULL
  # if(knots.option[2] == "quantile") z.knots <- select.knots(z,knots = k[2] -3)
  
  # object <- smooth.construct(z.object,
          # data = data.frame(z = data$z),knots=list(z=z.knots))
 
  z.knots <- construct.knots(data$z,knots[2],knots.option[2],p[2])
  z.List <- pspline.setting(data$z,knots=z.knots,p[2],m[2],
                            type="simple",knots.option=knots.option[2])
  
  Bz <- z.List$B
  Pz <- z.List$P
  cz <- nrow(Pz)
  
  z.object = list(z.List=z.List,
                  knots = z.knots,
                  knots.option = knots.option[2])
  
  B <- Matrix(kr(as.matrix(Bz),as.matrix(Bx)))
  P <- list()
  P[[1]] <- Matrix(kronecker(diag(cz),Px))
  P[[2]] <- Matrix(kronecker(Pz,diag(cx)))
  
  if(is.null(lambda)){
    fit <- igcv.wrapper(Y=y,X=B,P=P,N2=N2,
                       lower=lower,upper=upper,search.length=search.length)
    lambda <- fit$lambda
  }
  
  if(!is.null(lambda)) lambda <- rep(lambda,2)[1:2]
  
  Ptotal <- P[[1]]*lambda[1] + P[[2]]*lambda[2] 
  # print(lambda)
  temp <- tcrossprod(solve(crossprod(B) + Ptotal),B)
  theta <- temp%*%y
  Theta <- matrix(theta,nrow=ncol(Bx),ncol=ncol(Bz))
  
  
  res <- list(fitted.values = as.vector(B%*%theta), B = B,
              theta = theta,lambda = lambda,
              Theta=Theta,x.object = x.object,z.object=z.object, p = p)
  class(res) <- "bps"
  return(res)
}
