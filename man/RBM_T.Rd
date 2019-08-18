\name{RBM_T}
\alias{RBM_T}

\title{RBM_T: a R function for microarray and RNA-Seq data analysis for two-group comparisons
}
\description{ Use A Resampling-Based Empirical Bayes Approach to Assess Differential Expression or Identify differntially methylated loci in Two-Color Microarrays and RNA-Seq data sets.  
}
\usage{
RBM_T(aData, vec_trt, repetition, alpha)
}

\arguments{
  \item{aData}{
   The input data set with rows and columns denoting features and samples, respectively
}
  \item{vec_trt}{
   A vector for group notation such as 1s denote treatment group and 0s denote control group
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
RBM_T produces a named list with the following components:
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
The \code{\link{RBM_F}} function defined in this package.
The limma and marray packages.
}

\examples{
normal_data <- matrix(rnorm(200*6), 200, 6)
mydesign <- c(0,0,0,1,1,1)
norm_result <- RBM_T(normal_data,mydesign,50,0.05)
    
unif_data <- matrix(runif(200*7, 0.10, 0.95), 200, 7)
mydesign2 <- c(0,0,0, 1,1,1,1)
unif_result <- RBM_T(unif_data,mydesign2,100,0.05)
}

\keyword{ Resampling }
\keyword{ Microarray and RNA-Seq }
