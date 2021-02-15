\name{mface.sparse}
\alias{mface.sparse}

\title{
Fast covariance estimation for multivariate sparse functional data
}
\description{
The function is to estimate the mean and covariance function
from a cluster of multivariate functions/longitudinal observations.
}
\usage{
mface.sparse(data, newdata = NULL, 
             center = TRUE, argvals.new = NULL, knots = 7, 
             knots.option = "equally-spaced", 
             p = 3, m = 2, 
             lambda = NULL, lambda_mean = NULL, lambda_bps = NULL, 
             search.length = 14, lower = -3, upper = 10, 
             calculate.scores = FALSE, pve = 0.99)
}

\arguments{
  \item{data}{
  a list containing all functional outcomes.
  Each element is a data frame with three arguments:
  (1) \code{argvals}: observation times;
  (2) \code{subj}: subject indices;
  (3) \code{y}: values of observations for each dimension.
  Missing values not allowed.
}
  \item{newdata}{
  of the same strucutre as \code{data}; defaults to NULL, then no prediction. 
}
  \item{center}{
  logical. If TRUE, then Pspline smoothing of the population mean will be conducted and subtracted from the data before covariance smoothing; if FALSE, then the population mean will be just 0s.
}
  \item{argvals.new}{
  a vector of observation time points to evaluate mean function, covariance function, error variance and etc. .
}
  \item{knots}{
 the number of knots for B-spline basis functions to be used; defaults to 7.
    The resulting number of basis functions is the number of interior knots plus the degree of B-splines.
}
  \item{knots.option}{
  if \code{knots} specifies the number of knots, then \code{knots.option} will be used.  Default "equally-spaced", then equally-spaced knots in the range of observed time points will be selected; alternatively, "quantile": quantiles of the observed time points will be selected; see details.
}
  \item{p}{
  the degrees of B-splines; defaults to 3.
}
  \item{m}{
  the order of differencing penalty; defaults to 2.
}
  \item{lambda}{
  the value of the smoothing parameter for auto-covariance smoothing; defaults to NULL.
}
  \item{lambda_mean}{
  the value of the smoothing parameter for mean smoothing; defaults to NULL.
}
  \item{lambda_bps}{
  the value of the smoothing parameter for cross-covariance smoothing; defaults to NULL.
}
  \item{search.length}{
the number of equidistant (log scale) smoothing parameters to search; defaults to 14.
}

\item{lower, upper}{
bounds for log smoothing parameter; defaults are -3 and 10, respectively.
}
  \item{calculate.scores}{
if TRUE, scores will be calculated.
}
  \item{pve}{
Defaults 0.99. To select the number of eigenvalues by percentage of variance.
}
}
\details{
This is a generalized version of bivariate P-splines (Eilers and Marx, 2003) for
covariance smoothing of multivariate sparse functional or longitudinal data.
It uses tensor product B-spline basis functions and
employs differencing penalties on the assosciated parameter matrix.
The smoothing parameters in the method are selected
by leave-one-subject-out cross validation and is implemented with a fast algorithm.

If \code{center} is TRUE, then the population means will be calculated and are smoothed by
univariate P-spline smoothing: \code{pspline} (Eilers and Marx, 1996). This univariate
smoothing uses leave-one-subject-out cross validation to select the smoothing parameter.

If knots.option is "equally-spaced", then the differencing penalty in Eilers and Marx (2003) is used; if knots.option is "quantile" then the integrated squared second order derivative penalty in Wood (2016) is used.
                                                
}
\value{
\item{fit}{Univariate FPCA fit for each function}
\item{y.pred,mu.pred,Chat.diag.pred,se.pred,var.error.pred}{Predicted/estimated objects at \code{newdata$argvals}}
\item{Theta}{Estimated parameter matrix}
\item{argvals.new}{Vector of time points to evaluate population parameters}
\item{Chat.new,Cor.new,Cor.raw.new,Chat.raw.diag.new,var.error.new}{Estimated objects at \code{argvals.new}}
\item{eigenfunctions, eigenvalues}{Estimated eigenfunctions (scaled eigenvector) and eigenvalues at \code{argvals.new}}
\item{var.error.hat}{Estimated objects for each outcome}
\item{calculate.scores,scores}{if \code{calculate.scores} is TRUE (default to FALSE) and \code{newdata} is not NULL, then predicted scores will be calculated.}
\item{\dots}{\dots}
}
\references{
Cai Li, Luo Xiao, and Sheng Luo, 2020. Fast covariance estimation for multivariate sparse functional data. Stat, 9(1), p.e245, 
\doi{10.1002/sta4.245}.

Luo Xiao, Cai Li, William Checkley and Ciprian Crainiceanu, Fast covariance estimation
for sparse functional data, Stat. Comput., \doi{doi:10.1007/s11222-017-9744-8}.

Paul Eilers and Brian  Marx, Multivariate calibration with temperature
interaction using two-dimensional penalized signal regression,
Chemometrics and Intelligent Laboratory Systems 66 (2003), 159-174.

Paul Eilers and Brian Marx, Flexible smoothing with B-splines and penalties,
Statist. Sci., 11, 89-121, 1996.
 
Simon N. Wood, P-splines with derivative based penalties and tensor product 
  smoothing of unevenly distributed data, Stat. Comput., \doi{10.1007/s11222-016-9666-x}.
}
\author{
Cai Li <cli9@ncsu.edu> and Luo Xiao <lxiao5@ncsu.edu>
}

\seealso{
\code{face.sparse} in \code{face}
}

\examples{
\dontrun{
require(fields)

## settings
n <- 400
M <- 5
SNR <- 5
seed <- 118

set.seed(seed)

## data generation
if(M==5)  {
  N1 <- sample(3:7,n,replace=TRUE)
  N2 <- sample(3:7,n,replace=TRUE)
  N3 <- sample(3:7,n,replace=TRUE)
}

subj1 <- c()
subj2 <- c()
subj3 <- c()
for(i in 1:n){
  subj1 <- c(subj1,rep(i, N1[i]))
  subj2 <- c(subj2,rep(i, N2[i]))
  subj3 <- c(subj3,rep(i, N3[i]))
}
t1 <- runif(sum(N1))
t2 <- runif(sum(N2))
t3 <- runif(sum(N3))
tnew <- seq(0,1,length=100)
y1 <- 5*sin(2*pi*t1)
y2 <- 5*cos(2*pi*t2)
y3 <- 5*(t3-1)^2

x1 <- t(matrix(rep(5*sin(2*pi*tnew),n),length(tnew),n)) 
x2 <- t(matrix(rep(5*cos(2*pi*tnew),n),length(tnew),n)) 
x3 <- t(matrix(rep(5*(tnew-1)^2,n),length(tnew),n)) 

psi11 <- function(x){sqrt(2/3)*sin(2*pi*x)}
psi12 <- function(x){sqrt(2/3)*cos(4*pi*x)}
psi13 <- function(x){sqrt(2/3)*sin(4*pi*x)}

psi21 <- function(x){sqrt(2/3)*sin((1-1/2)*pi*x)}
psi22 <- function(x){sqrt(2/3)*sin((2-1/2)*pi*x)}
psi23 <- function(x){sqrt(2/3)*sin((3-1/2)*pi*x)}

psi31 <- function(x){sqrt(2/3)*sin(1*pi*x)}
psi32 <- function(x){sqrt(2/3)*sin(2*pi*x)}
psi33 <- function(x){sqrt(2/3)*sin(3*pi*x)}

Lambda <- c(2,1,0.5)*3

x <- matrix(NA,nrow=n*length(tnew),ncol=3)
xi <- matrix(NA,nrow=n,ncol=3)
for(k in 1:3){xi[,k] = rnorm(n)*sqrt(Lambda[k])}

for(i in 1:n){
  seq1 <- (sum(N1[1:i])-N1[i]+1):(sum(N1[1:i]))
  seq2 <- (sum(N2[1:i])-N2[i]+1):(sum(N2[1:i]))
  seq3 <- (sum(N3[1:i])-N3[i]+1):(sum(N3[1:i]))
  
  Xt = xi[i,1]*c(psi11(t1[seq1]),psi21(t2[seq2]),psi31(t3[seq3])) + 
    xi[i,2]*c(psi12(t1[seq1]),psi22(t2[seq2]),psi32(t3[seq3])) + 
    xi[i,3]*c(psi13(t1[seq1]),psi23(t2[seq2]),psi33(t3[seq3]))
  
  y1[seq1] = y1[seq1] + Xt[1:N1[i]]
  y2[seq2] = y2[seq2] + Xt[N1[i]+1:N2[i]]
  y3[seq3] = y3[seq3] + Xt[N1[i]+N2[i]+1:N3[i]]
  
  x[((i-1)*length(tnew)+1) :(length(tnew)*i),] = c(x1[i,], x2[i,], x3[i,]) + 
    xi[i,1]*c(psi11(tnew),psi21(tnew),psi31(tnew)) + 
    xi[i,2]*c(psi12(tnew),psi22(tnew),psi32(tnew)) + 
    xi[i,3]*c(psi13(tnew),psi23(tnew),psi33(tnew))
}

True_C <- Lambda[1]*c(psi11(tnew),psi21(tnew),psi31(tnew))\%x\%
  t(c(psi11(tnew), psi21(tnew), psi31(tnew))) +  
  Lambda[2]*c(psi12(tnew), psi22(tnew), psi32(tnew))\%x\%
  t(c(psi12(tnew), psi22(tnew), psi32(tnew))) + 
  Lambda[3]*c(psi13(tnew), psi23(tnew), psi33(tnew))\%x\%
  t(c(psi13(tnew), psi23(tnew), psi33(tnew)))

sigma <- sqrt(sum(Lambda)/(3*SNR))

## observed data
y1 <- y1 + rnorm(sum(N1))*sigma
y2 <- y2 + rnorm(sum(N2))*sigma
y3 <- y3 + rnorm(sum(N3))*sigma

# true trajectories
x1 <- t(matrix(x[,1],length(tnew),n))
x2 <- t(matrix(x[,2],length(tnew),n))
x3 <- t(matrix(x[,3],length(tnew),n))

true_eigenfunctions <- eigen(True_C)$vectors*sqrt(length(tnew))
true_eigenvalues <- eigen(True_C)$values/length(tnew)

## organize data and apply mFACEs
data <- list("y1" = data.frame("subj"= subj1, "argvals" = t1, "y" = y1),
             "y2" = data.frame("subj"= subj2, "argvals" = t2, "y" = y2),
             "y3" = data.frame("subj"= subj3, "argvals" = t3, "y" = y3)) 
fit <- mface.sparse(data, argvals.new = tnew)

## heat maps of correlation function
par(mfrow=c(1,2),mar=c(5,2.5,3,5.5))
var <- diag(True_C)
var[var==0] = 1e-6
diag(True_C) = var
cor <- diag(1/sqrt(var)) \%*\% True_C \%*\% diag(1/sqrt(var))
image.plot(cor,
           xlab="t",
           ylab="s",main=bquote("Truth"),
           cex.axis=2.0,cex.lab=2.0,cex.main=2.0,legend.cex=2.0,axes=F,
           zlim=c(-1,1))
axis(1, at=seq(0,1,length.out = 10),labels=c(rep(c("0","0.25","0.75"),3),"1"))
axis(2, at=seq(0,1,length.out = 10),labels=c(rep(c("0","0.25","0.75"),3),"1"))

cor <- as.matrix(fit$Cor.new)
image.plot(cor,
           xlab="t",
           ylab="s",main=bquote("Estimated"),
           cex.axis=2.0,cex.lab=2.0,cex.main=2.0,legend.cex=2.0,axes=F,
           zlim=c(-1,1))
axis(1, at=seq(0,1,length.out = 10),labels=c(rep(c("0","0.25","0.75"),3),"1"))
axis(2, at=seq(0,1,length.out = 10),labels=c(rep(c("0","0.25","0.75"),3),"1"))

## prediction of several subjects
for(i in 1:3){
  sel <- lapply(data, function(x){which(x$subj==i)})
  dat_i <- mapply(function(data, sel){data[sel,]}, 
                  data = data, sel = sel, SIMPLIFY = F)
  dat_i_pred <- lapply(dat_i, function(x){
    data.frame(subj=rep(x$subj[1],nrow(x) + length(tnew)),
               argvals = c(rep(NA,nrow(x)),tnew),
               y = rep(NA,nrow(x) + length(tnew)))
  })
  for(j in 1:length(dat_i)){
    dat_i_pred[[j]][1:nrow(dat_i[[j]]), ] <- dat_i[[j]]
  }
  pred <- predict(fit, dat_i_pred)
  y_pred <- mapply(function(pred_y.pred, dat_i){
    pred_y.pred[nrow(dat_i)+1:length(tnew)]}, pred_y.pred = pred$y.pred, 
    dat_i = dat_i, SIMPLIFY = T)

  pre <- pred
  par(mfrow=c(1,3),mar=c(5,5.8,3,1))
  Ylim = c(-12,12)
  Xlim = c(0,1)
  Ylab = bquote(y^(1))
  Xlab = "t"
  main = paste("Subject ", dat_i[[1]][1,1],sep="")
  idx = (nrow(dat_i[[1]])+1):(nrow(dat_i[[1]])+length(tnew))
  plot(dat_i[[1]][,"argvals"],dat_i[[1]][,"y"],ylim=Ylim,xlim=Xlim,ylab=Ylab,xlab=Xlab,
       main=main,cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0,pch=1)
  lines(tnew,pre$y.pred$y1[idx],col="red",lwd=2)
  lines(tnew,pre$y.pred$y1[idx]-1.96*pre$se.pred$y1[idx],col="blue",lwd=2,lty=2)
  lines(tnew,pre$y.pred$y1[idx]+1.96*pre$se.pred$y1[idx],col="blue",lwd=2,lty=2)
  lines(tnew,x1[i,],col="purple",lwd=2)

  Ylab = bquote(y^(2))
  Xlab = "t"
  main = paste("Subject ", dat_i[[1]][1,1],sep="")
  idx = (nrow(dat_i[[2]])+1):(nrow(dat_i[[2]])+length(tnew))
  plot(dat_i[[2]][,"argvals"],dat_i[[2]][,"y"],ylim=Ylim,xlim=Xlim,ylab=Ylab,xlab=Xlab,
       main=main,cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0,pch=1)
  lines(tnew,pre$y.pred$y2[idx],col="red",lwd=2)
  lines(tnew,pre$y.pred$y2[idx]-1.96*pre$se.pred$y2[idx],col="blue",lwd=2,lty=2)
  lines(tnew,pre$y.pred$y2[idx]+1.96*pre$se.pred$y2[idx],col="blue",lwd=2,lty=2)
  lines(tnew,x2[i,],col="purple",lwd=2)

  Ylab = bquote(y^(3))
  Xlab = "t"
  main = paste("Subject ", dat_i[[1]][1,1],sep="")
  idx = (nrow(dat_i[[3]])+1):(nrow(dat_i[[3]])+length(tnew))
  plot(dat_i[[3]][,"argvals"],dat_i[[3]][,"y"],ylim=Ylim,xlim=Xlim,ylab=Ylab,xlab=Xlab,
       main=main,cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0,pch=1)
  lines(tnew,pre$y.pred$y3[idx],col="red",lwd=2)
  lines(tnew,pre$y.pred$y3[idx]-1.96*pre$se.pred$y3[idx],col="blue",lwd=2,lty=2)
  lines(tnew,pre$y.pred$y3[idx]+1.96*pre$se.pred$y3[idx],col="blue",lwd=2,lty=2)
  lines(tnew,x3[i,],col="purple",lwd=2)

}

}
}


\keyword{~mface.sparse}
