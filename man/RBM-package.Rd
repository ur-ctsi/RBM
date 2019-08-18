\name{RBM-package}
\alias{RBM-package}
\alias{RBM}
\docType{package}
\title{RBM:a package for microarray and RNA-Seq data analysis}
\description{
Use A Resampling-Based Empirical Bayes Approach to Assesse Differential Expression or Identifying differntially methylated loci in Two-Color Microarrays and RNA-Seq data sets. Significant features selected through RBM_T or RBM_F functions could be further used as input for pathway analysis or experimental vilidations.
}

\details{
\tabular{ll}{
Package: \tab RBM\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2014-10-05\cr
Depends: \tab R (>= 3.0.0), limma, marray \cr
License: \tab GPL (>= 2) \cr
}
}

\author{
Dongmei Li and Chin-Yuan Liang
Maintainer: Dongmei Li <dongmeiliur@gmail.com> and Chin-Yuan Liang <liang.tony@gmail.com>
}

\references{
Li D, Le Pape MA, Parikh NI, Chen WX, Dye TD (2013) Assessing Differential Expression in Two-Color Microarrays: A Resampling-Based Empirical Bayes Approach. PLoS ONE 8(11): e80099. doi: 10.1371/journal.pone.0080099
}

\keyword{Resampling, Empirical Bayes, Microarray, RNA-Seq}

\seealso{
The \code{\link{RBM_T}} and \code{\link{RBM_F}} functions defined in this package.
The limma and marray packages.
}

\examples{
normal_data <- matrix(rnorm(200*6), 200, 6)
mydesign <- c(0,0,0,1,1,1)
norm_result <- RBM_T(normal_data,mydesign,50,0.05)
    
unif_data <- matrix(runif(200*7, 0.10, 0.95), 200, 7)
mydesign2 <- c(0,0,0, 1,1,1,1)
unif_result <- RBM_T(unif_data,mydesign2,100,0.05)
    
normdata_F <- matrix(rnorm(200*9, 0, 2), 200, 9)   
mydesign_F <- c(0, 0, 0, 1, 1, 1, 2, 2, 2)
aContrast <- c("X1-X0", "X2-X1", "X2-X0")
normresult_F <- RBM_F(normdata_F, mydesign_F, aContrast, 100, 0.05) 
     
unifdata_F <- matrix(runif(200*18, 0.15, 0.98), 200, 18) 
mydesign2_F <- c(rep(0, 6), rep(1, 6), rep(2, 6))
aContrast <- c("X1-X0", "X2-X1", "X2-X0")
unifresult_F <- RBM_F(unifdata_F, mydesign2_F, aContrast, 100, 0.05) 
}


