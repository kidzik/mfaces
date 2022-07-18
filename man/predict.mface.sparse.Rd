\name{predict.mface.sparse}
\alias{predict.mface.sparse}
\title{
Subject-specific curve prediction from a mface.sparse fit
}
\description{
Predict subject-specific curves based on a fit from "mface.sparse".
}
\usage{
\method{predict}{mface.sparse}(object, newdata, calculate.scores = T, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted object from the R function  "mface.sparse".
}
  \item{newdata}{
a list containing all functional outcomes. Each element is a data frame with three arguments: 
(1) \code{argvals}: observation times;
(2) \code{subj}: subject indices;
(3) \code{y}: values of observations for each dimension.
NA values are allowed in "y" but not in the other two. 
}
  \item{calculate.scores}{
if TRUE, scores will be calculated.
}

  \item{\dots}{
further arguments passed to or from other methods.
}
}

\details{
This function makes prediction based on observed data for each subject. So for each subject,
it requires at least one observed data. For the time points prediction is desired but no observation is available, just make the corresponding \code{data$y} as NA.
}
\value{
\item{object}{A "mface.sparse" fit}
\item{newdata}{Input data}
\item{y.pred,mu.pred,se.pred,Chat.diag.pred,var.error.pred}{Predicted/estimated objects at the observation time points in \code{newdata}}
\item{rand_eff}{if \code{calculate.scores} in \code{object} is TRUE (typically TRUE), then predicted scores \code{rand_eff$scores} will be calculated.}
\item{\dots}{\dots}
}

\references{
Cai Li, Luo Xiao, and Sheng Luo, 2020. Fast covariance estimation for multivariate sparse functional data. Stat, 9(1), p.e245, 
\doi{10.1002/sta4.245}.
}
\author{
Cai Li <cli9@ncsu.edu>
}


\examples{
# See the examples for "mface.sparse".
}

\keyword{~mface.sparse}
\keyword{~prediction}

