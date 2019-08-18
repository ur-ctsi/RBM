\name{RBM_F}
\alias{RBM_F}

\title{RBM_F: a R function for microarray and RNA-Seq data analysis for designs with more than two groups
}
\description{Use A Resampling-Based Empirical Bayes Approach to Assess Differential Expression in Two-Color Microarrays and RNA-Seq data sets for designs with more than two groups. 
}
\usage{
RBM_F(aData, vec_trt, aContrast, repetition, alpha)
}

\arguments{
  \item{aData}{
    The input data set with rows and columns denoting features and samples, respectively
  }
  \item{vec_trt}{
    A vector for group notation such as 1s denote treatment group and 0s denote control group
  }
  \item{aContrast}{
   A vector for contrast. For example: if we want to compare group 1 with group 0, group 2 with group 1, and group 2 with group 0, then the contrast vector will be ("X1-X0", "X2"-"X1", "X2-X0")
  }
  \item{repetition}{
   The number of resamplings used in the analysis. You could use 1000 or higher number
  }
  \item{alpha}{
   The signifiance level
  }
}

\details{Combine resampling with empirical Bayes approach for Microarrays and RNA-Seq data analysis.}

\value{
RBM_F produces a named list with the following components:
  \item{ordfit_t}{orignal t statistics}             
  \item{ordfit_pvalue}{original p-values from lmFit and eBayes}        
  \item{ordfit_beta0}{estimated mean for the control group}        
  \item{ordfit_beta1}{estimated mean difference between treatment and control group}         
  \item{permutation_p}{calculated p-values from permutation method based on resampled test statistics}  
  \item{bootstrap_p}{calculated p-values from bootstrap method based on resampled test statistics}
}

\references{
Li D, Le Pape MA, Parikh NI, Chen WX, Dye TD (2013) Assessing Differential Expression in Two-Color Microarrays: A Resampling-Based Empirical Bayes Approach. PLoS ONE 8(11): e80099. doi: 10.1371/journal.pone.0080099
}
\author{Dongmei Li and Chin-Yuan Liang}

\seealso{
The \code{\link{RBM_T}} function defined in this package.
The limma and marray packages.
}

\examples{
normdata_F <- matrix(rnorm(200*9, 0, 2), 200, 9)   
mydesign_new <- c(0, 0, 0, 1, 1, 1, 2, 2, 2)
aContrast <- c("X1-X0", "X2-X1", "X2-X0")
normresult_F <- RBM_F(normdata_F, mydesign_new, aContrast, 100, 0.05) 
     
unifdata_F <- matrix(runif(200*18, 0.15, 0.98), 200, 18) 
mydesign2_new <- c(rep(0, 6), rep(1, 6), rep(2, 6))
aContrast <- c("X1-X0", "X2-X1", "X2-X0")
unifresult_F <- RBM_F(unifdata_F, mydesign2_new, aContrast, 100, 0.05) 
}

\keyword{Resampling }
\keyword{Microarray and RNA-Seq }
