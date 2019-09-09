## Set the working directory
setwd("\\path_to_your_working_dir\\")

## Load the functions needed for simulation analysis
source("Simulation_functions.R")

### This script is for simulating scenarios where there are 1, 2, 3, or 4 regulatory SNPs.
### We calculate "predicted" R2 values using the formulas listed in Figure 1 (1 or 2 regulatory SNPs) or Table 1 (3 or 4 regulatory SNPs) in the main text.
### Lastly we produce scatter plots to visualize correlations between predicted R2 and estimated R2 values. "Estimated" R2 values = R2 values determined from 
### linear regression analysis of simulated mRNA expression levels vs. simulated SNP genotype data

#########################################################
############## One regulatory SNP simulation ############
#########################################################

vb <- 0.8  # Variance explained by the regulatory SNP
N <- 3000  # Number of simulations
res <- data.frame(R2A = rep(0, N), R2A.m = 0, pa= 0, pb=0, rAB=0)

# For each simulation, calculate the estimated R2 and predicted R2
l <- 1
for (i in 1:N) {
  res[l, ] <- simu_sf(vb, err = 0.2)  ## error ~ Norm(0, 0.2 * total_variance)
  l <- l + 1
}
res <- res[!is.na(res$R2A.m) & (res$R2A.m != Inf), ]
summary(lm(res[, 1]~res[, 2]))$r.squared

## Plot predicted R2 vs. estimated R2
pdf("figureS1-1rsnp-with-err.pdf", height=6, width=6)
plot(res[, 2], res[, 1], xlab = NA, ylab = NA, cex.axis = 2, col = "skyblue")
dev.off()

#########################################################
############## Two regulatory SNP simulation ############
#########################################################
vb <- 0.5
vc <- 0.4
N <- 3000
res <- data.frame(R2A = rep(0, N), R2A.m = 0, 
                  pa= 0, pb=0, pc=0, rAB=0, rBC=0, rAC=0, rBP=0, rCP=0,
                  v=0, covbc=0)
l <- 1
for (i in 1:N) {
  res[l, ] <- simu2_sf(vb, vc, err = 0.2) ## error ~ Norm(0, 0.2 * total_variance)
  l <- l + 1
}
res <- res[!is.na(res$R2A.m) & (res$R2A.m != Inf), ]
summary(lm(res[, 1]~res[, 2]))$r.squared

## Plot predicted R2 vs. estimated R2
pdf("figureS1-2rsnp-with-err.pdf", height=6, width=6)
plot(res[, 2], res[, 1], xlab = NA, ylab = NA, cex.axis = 2, col = "skyblue")
dev.off()

###########################################################
############## Three regulatory SNP simulation ############
###########################################################
vb <- 0.4
vc <- 0.3
vd <- 0.2
N <- 3000
res <- data.frame(R2A = rep(0, N), R2A.m = 0, pa = 0, pb = 0, pc = 0, pd = 0, rAB = 0, rBC = 0, 
                  rAC = 0, rAD = 0, rBD = 0, rCD = 0, rBP = 0, rCP = 0, rDP = 0)
l <- 1
for (i in 1:N) {
  res[l, ] <- simu3_sf(vb, vc, vd, err = 0.2) ## error ~ Norm(0, 0.2 * total_variance)
  l <- l + 1
}
res <- res[!is.na(res$R2A.m) & (res$R2A.m != Inf), ]
summary(lm(res[, 1]~res[, 2]))$r.squared

## Plot predicted R2 vs. estimated R2
pdf("figureS1-3rsnp-with-err.pdf", height=6, width=6)
plot(res[, 2], res[, 1], xlab = NA, ylab = NA, cex.axis = 2, col = "skyblue")
dev.off()

##########################################################
############## Four regulatory SNP simulation ############
##########################################################

vb <- 0.4
vc <- 0.6
vd <- 0.2
ve <- 0.3
N <- 3000
res <- data.frame(R2A = rep(0, N), R2A.m = 0) ## error ~ Norm(0, 0.2 * total_variance)
l <- 1
for (i in 1:N) {
  if (i %% 100 == 0) {print(i)}
  res[l, ] <- simu2_sf(vb, vc, err=0.2)
  l <- l + 1
}
res <- res[!is.na(res$R2A.m) & (res$R2A.m != Inf), ]
summary(lm(res[, 1]~res[, 2]))$r.squared

## Plot predicted R2 vs. estimated R2
pdf("figureS1-4rsnp-with-err.pdf", height=6, width=6)
plot(res[, 2], res[, 1], xlab = NA, ylab = NA, cex.axis = 2, col = "skyblue")
dev.off()
##########################################################################
