RBM_T<-
  function (aData, vec_trt, repetition, alpha) 
  {
    mydata <- as.matrix(aData)
    rowsdim <- dim(aData)[1]
    colsdim <- dim(aData)[2]
    mydesign <- cbind(1, vec_trt)
    colnames(mydesign) <- c("Grp1", "Grp2vs1")
    tmp <- interRBM(aData, , mydesign)
    ordfit_t <- tmp$tstat
    ordfit_pvalue <- tmp$pvalue
    ordfit_beta0 <- tmp$beta0
    ordfit_beta1 <- tmp$beta1
    permutation_t <- matrix(NA, rowsdim, repetition)
    permutation_pvalue <- matrix(NA, rowsdim, repetition)
    for (t in 1:repetition) {
      aPermutation <- sample(1:colsdim, replace = FALSE)
      tmp <- interRBM(aData, aPermutation, mydesign)
      permutation_t[, t] <- tmp$tstat
      permutation_pvalue[, t] <- tmp$pvalue
    }
    bootstrap_t <- matrix(NA, rowsdim, repetition)
    bootstrap_pvalue <- matrix(NA, rowsdim, repetition)
    for (t in 1:repetition) {
      aBootstrap <- sample(1:colsdim, replace = TRUE)
      tmp <- interRBM(aData, aBootstrap, mydesign)
      bootstrap_t[, t] <- tmp$tstat
      bootstrap_pvalue[, t] <- tmp$pvalue
    }
    tmp <- interRBM_eval(ordfit_t, ordfit_pvalue, permutation_t, 
                         permutation_pvalue, repetition)
    permutation_p_tstat <- tmp$p_tstat
    permutation_p_pvalue <- tmp$p_pvalue
    tmp <- interRBM_eval(ordfit_t, ordfit_pvalue, bootstrap_t, 
                         bootstrap_pvalue, repetition)
    bootstrap_p_tstat <- tmp$p_tstat
    bootstrap_p_pvalue <- tmp$p_pvalue
    mylist <- list(ordfit_t = ordfit_t, ordfit_pvalue = ordfit_pvalue, 
                   ordfit_beta0 = ordfit_beta0, ordfit_beta1 = ordfit_beta1, 
                   permutation_p = permutation_p_tstat, bootstrap_p = bootstrap_p_tstat)
    return(mylist)
  }

interRBM<-
  function (aData, aList, aDesign) 
  {
    if (exists("aList")) {
      mydata <- aData[, aList]
    }
    else {
      mydata <- aData
    }
    myfit <- lmFit(mydata, aDesign)
    myfit <- eBayes(myfit)
    tstat <- myfit$t[, 2]
    pvalue <- myfit$p.value[, 2]
    beta0 <- myfit$coefficients[, 1]
    beta1 <- myfit$coefficients[, 2]
    mylist <- list(tstat = tstat, pvalue = pvalue, beta0 = beta0, 
                   beta1 = beta1)
    return(mylist)
  }

RBM_F<-function(aData, vec_trt, aContrast, repetition, alpha)
{
  mydata<-as.matrix(aData)
  rowsdim<-dim(aData)[1]
  colsdim<-dim(aData)[2]
  mydesign_F<-model.matrix(~0+factor(vec_trt))
  colnames(mydesign_F)<-make.names(levels(factor(vec_trt)))
  mycontrast<-aContrast
  tmp <- interRBM_F(aData, , mydesign_F, mycontrast)
  ordfit_t <- tmp$tstat
  ordfit_pvalue <- tmp$pvalue
  ordfit_beta1 <- tmp$beta1
  permutation_p_tstat<-permutation_p_pvalue<-bootstrap_p_tstat<-bootstrap_p_pvalue<-matrix(NA, nrow=rowsdim, ncol=length(aContrast))
  permutation_t <- matrix(NA, rowsdim, repetition)
  permutation_pvalue <- matrix(NA, rowsdim, repetition)
  bootstrap_t <- matrix(NA, rowsdim, repetition)
  bootstrap_pvalue <- matrix(NA, rowsdim, repetition)
  for (i in 1:length(aContrast))
  {
    for (t in 1:repetition) {
      aPermutation <- sample(1:colsdim, replace = FALSE)
      tmp <- interRBM_F(aData, aPermutation, mydesign_F, mycontrast)
      permutation_t[, t] <- tmp$tstat[, i]
      permutation_pvalue[, t] <- tmp$pvalue[, i]
    }
    
    for (t in 1:repetition) {
      aBootstrap <- sample(1:colsdim, replace = TRUE)
      tmp <- interRBM_F(aData, aBootstrap, mydesign_F, mycontrast)
      bootstrap_t[, t] <- tmp$tstat[, i]
      bootstrap_pvalue[, t] <- tmp$pvalue[, i]
    }
    tmp <- interRBM_eval(ordfit_t[, i], ordfit_pvalue[, i], permutation_t, 
                         permutation_pvalue, repetition)
    permutation_p_tstat[, i] <- tmp$p_tstat
    permutation_p_pvalue[, i] <- tmp$p_pvalue
    tmp <- interRBM_eval(ordfit_t[, i], ordfit_pvalue[, i], bootstrap_t, 
                         bootstrap_pvalue, repetition)
    bootstrap_p_tstat[, i] <- tmp$p_tstat
    bootstrap_p_pvalue[, i] <- tmp$p_pvalue
  }
  mylist <- list(ordfit_t = ordfit_t, ordfit_pvalue = ordfit_pvalue, ordfit_beta1 = ordfit_beta1, 
                 permutation_p = permutation_p_tstat, bootstrap_p = bootstrap_p_tstat)
  return(mylist)
}

interRBM_eval<-
  function (ord_vec_t, ord_vec_pvalue, eval_matrix_t, eval_matrix_pvalue, repetition) 
  {
    rowsdim <- length(ord_vec_t)
    p_pvalue <- rep(NA, rowsdim)
    p_tstat <- rep(NA, rowsdim)
    for (g in 1:rowsdim) {
      p_pvalue[g] <- sum(eval_matrix_pvalue[g, ] <= ord_vec_pvalue[g])/repetition
      p_tstat[g] <- sum(abs(eval_matrix_t[g, ]) >= abs(ord_vec_t[g]))/repetition
    }
    mylist <- list(p_pvalue = p_pvalue, p_tstat = p_tstat)
    return(mylist)
  }


interRBM_F<-function (aData, aList, aDesign, aContrast) 
{
  if (exists("aList")) {
    mydata <- aData[, aList]
  }
  else {
    mydata <- aData
  }
  myfit <- lmFit(mydata, aDesign)
  fpvalue<-fbeta1<-ftstat<-matrix(NA, nrow=dim(mydata)[1], ncol=length(aContrast))
  for (i in 1:length(aContrast))
  {
    contrast.matrix<-makeContrasts(aContrast[i], levels=aDesign)
    myfit2<-contrasts.fit(myfit, contrast.matrix)
    myfit2<-eBayes(myfit2)
    tstat <- myfit2$t
    pvalue <- myfit2$p.value
    beta1 <- myfit2$coefficients
    ftstat[, i]<-tstat
    fpvalue[, i]<-pvalue
    fbeta1[, i]<-beta1      
  }
  mylist <- list(tstat = ftstat, pvalue = fpvalue, beta1 = fbeta1)
  return(mylist)
}

