### R code from vignette source 'RBM.Rnw'

###################################################
### code chunk number 1: RBM.Rnw:47-50 (eval = FALSE)
###################################################
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install("RBM")


###################################################
### code chunk number 2: RBM.Rnw:56-57 (eval = FALSE)
###################################################
## library(RBM)


###################################################
### code chunk number 3: RBM.Rnw:70-86
###################################################
library(RBM)
normdata <- matrix(rnorm(1000*6, 0, 1),1000,6)
mydesign <- c(0,0,0,1,1,1)
myresult <- RBM_T(normdata,mydesign,100,0.05)
summary(myresult)
sum(myresult$permutation_p<=0.05)
  
which(myresult$permutation_p<=0.05)

sum(myresult$bootstrap_p<=0.05)

which(myresult$bootstrap_p<=0.05)
permutation_adjp <- p.adjust(myresult$permutation_p, "BH")
sum(permutation_adjp<=0.05)
bootstrap_adjp <- p.adjust(myresult$bootstrap_p, "BH")
sum(bootstrap_adjp<=0.05)


###################################################
### code chunk number 4: RBM.Rnw:90-98
###################################################
unifdata <- matrix(runif(1000*7,0.10, 0.95), 1000, 7)
mydesign2 <- c(0,0,0, 1,1,1,1)
myresult2 <- RBM_T(unifdata,mydesign2,100,0.05)
sum(myresult2$permutatioin_p<=0.05)
sum(myresult2$bootstrap_p<=0.05)
which(myresult2$bootstrap_p<=0.05)
bootstrap2_adjp <- p.adjust(myresult2$bootstrap_p, "BH")
sum(bootstrap2_adjp<=0.05)


###################################################
### code chunk number 5: RBM.Rnw:104-130
###################################################
normdata_F <- matrix(rnorm(1000*9,0,2), 1000, 9)
mydesign_F <- c(0, 0, 0, 1, 1, 1, 2, 2, 2)
aContrast <- c("X1-X0", "X2-X1", "X2-X0")

myresult_F <- RBM_F(normdata_F, mydesign_F, aContrast, 100, 0.05)
summary(myresult_F)
sum(myresult_F$permutation_p[, 1]<=0.05)
sum(myresult_F$permutation_p[, 2]<=0.05)
sum(myresult_F$permutation_p[, 3]<=0.05)

which(myresult_F$permutation_p[, 1]<=0.05)
which(myresult_F$permutation_p[, 2]<=0.05)
which(myresult_F$permutation_p[, 3]<=0.05)
    
con1_adjp <- p.adjust(myresult_F$permutation_p[, 1], "BH")
sum(con1_adjp<=0.05/3)
  
con2_adjp <- p.adjust(myresult_F$permutation_p[, 2], "BH")
sum(con2_adjp<=0.05/3)

con3_adjp <- p.adjust(myresult_F$permutation_p[, 3], "BH")
sum(con3_adjp<=0.05/3)

which(con2_adjp<=0.05/3)

which(con3_adjp<=0.05/3)


###################################################
### code chunk number 6: RBM.Rnw:134-161
###################################################
unifdata_F <- matrix(runif(1000*18, 0.15, 0.98), 1000, 18)
mydesign2_F <- c(rep(0, 6), rep(1, 6), rep(2, 6))
aContrast <- c("X1-X0", "X2-X1", "X2-X0")

myresult2_F <- RBM_F(unifdata_F, mydesign2_F, aContrast, 100, 0.05)
summary(myresult2_F)

sum(myresult2_F$bootstrap_p[, 1]<=0.05)
  
sum(myresult2_F$bootstrap_p[, 2]<=0.05)

sum(myresult2_F$bootstrap_p[, 3]<=0.05)

which(myresult2_F$bootstrap_p[, 1]<=0.05)

which(myresult2_F$bootstrap_p[, 2]<=0.05)
 
which(myresult2_F$bootstrap_p[, 3]<=0.05)
 
con21_adjp <- p.adjust(myresult2_F$bootstrap_p[, 1], "BH")
sum(con21_adjp<=0.05/3)

con22_adjp <- p.adjust(myresult2_F$bootstrap_p[, 2], "BH")
sum(con22_adjp<=0.05/3)

con23_adjp <- p.adjust(myresult2_F$bootstrap_p[, 3], "BH")
sum(con23_adjp<=0.05/3)


###################################################
### code chunk number 7: RBM.Rnw:170-194
###################################################
system.file("data", package = "RBM")
data(ovarian_cancer_methylation)
summary(ovarian_cancer_methylation)
ovarian_cancer_data <- ovarian_cancer_methylation[, -1]
label <- c(1, 1, 0, 0, 1, 1, 0, 0)
diff_results <- RBM_T(aData=ovarian_cancer_data, vec_trt=label, repetition=100, alpha=0.05)
summary(diff_results)
sum(diff_results$ordfit_pvalue<=0.05)
sum(diff_results$permutation_p<=0.05)
sum(diff_results$bootstrap_p<=0.05)

ordfit_adjp <- p.adjust(diff_results$ordfit_pvalue, "BH")
sum(ordfit_adjp<=0.05)
perm_adjp <- p.adjust(diff_results$permutation_p, "BH")
sum(perm_adjp<=0.05)
boot_adjp <- p.adjust(diff_results$bootstrap_p, "BH")
sum(boot_adjp<=0.05)

diff_list_perm <- which(perm_adjp<=0.05)
diff_list_boot <- which(boot_adjp<=0.05)
sig_results_perm <- cbind(ovarian_cancer_methylation[diff_list_perm, ], diff_results$ordfit_t[diff_list_perm], diff_results$permutation_p[diff_list_perm])
print(sig_results_perm)
sig_results_boot <- cbind(ovarian_cancer_methylation[diff_list_boot, ], diff_results$ordfit_t[diff_list_boot], diff_results$bootstrap_p[diff_list_boot])
print(sig_results_boot)


