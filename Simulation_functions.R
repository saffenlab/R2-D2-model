
## Simulate one regulatory SNPs
simu_sf <- function (vb, err = 0) {

  p.haplo <- random_4()
  pAB <- p.haplo[1]
  pAb <- p.haplo[2]
  paB <- p.haplo[3]
  pab <- p.haplo[4]
  
  pa <- sum(p.haplo[as.logical(c(0, 0, 1, 1))])
  pb <- sum(p.haplo[as.logical(c(0, 1, 0, 1))])
  
  pDiplo <- c(pAB*pAB, pAB*pAb, pAb*pAB, pAb*pAb,
              pAB*paB, pAB*pab, pAb*paB, pAb*pab, paB*pAB, paB*pAb, pab*pAB, pab*pAb, 
              paB*paB, paB*pab, pab*paB, pab*pab)
  
  N <- 1000 # sample size
  nDiplo <- round(N*pDiplo)
  nind <- sum(nDiplo)
  
  codeA <- c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2)
  codeB <- c(0, 1, 1, 2, 0, 1, 1, 2, 0, 1, 1, 2, 0, 1, 1, 2)
  
  genoA <- rep(0, nind)
  genoB <- rep(0, nind)
  
  l <- 1
  for (i in 1:16) {
    if (nDiplo[i]) {
      genoA[l:(l+nDiplo[i]-1)] <- rep(codeA[i], nDiplo[i])
      genoB[l:(l+nDiplo[i]-1)] <- rep(codeB[i], nDiplo[i])
      l <- l + nDiplo[i]
    }
  }
  
  gg.b <- eogCalc(pb, vb)
  if (!err)
    pheno <- sapply(1:nind, function(x) {gg.b[genoB[x]+1]})
  else {
    pheno <- sapply(1:nind, function(x) {gg.b[genoB[x]+1] + rnorm(1, 0, sqrt(vb)*err)})
  }
  
  R2A <- summary(lm(pheno~genoA))$r.squared
  R2B <- summary(lm(pheno~genoB))$r.squared
  
  rAB <- cor(genoA, genoB)
  R2A.m <- R2B*(rAB)^2
  
  return(c(R2A, R2A.m, pa, pb, rAB))
}

## Simulate two regulatory SNPs
simu2_sf <- function(vb, vc, err) {
  
  code.abc <- data.frame(a = c(0, 0, 1, 1, 0, 0, 1, 1), b = c(0, 1, 0, 1, 0, 1, 0, 1), 
                         c = c(0, 1, 1, 0, 1, 0, 0, 1))
  p.haplo <- as.numeric(random_8()) ## the only difference with simu()
  pa <- sum(p.haplo[as.logical(c(0, 0, 1, 1, 0, 0, 1, 1))])
  pb <- sum(p.haplo[as.logical(c(0, 1, 0, 1, 0, 1, 0, 1))])
  pc <- sum(p.haplo[as.logical(c(0, 1, 1, 0, 1, 0, 0, 1))])
  
  n.haplo <- length(p.haplo)
  p.diplo <- rep(0, n.haplo^2)
  n.diplo <- length(p.diplo)
  code.a <- rep(0, n.diplo)
  code.b <- rep(0, n.diplo)
  code.c <- rep(0, n.diplo)
  
  l <- 1
  for (i in 1:n.haplo) 
    for (j in 1:n.haplo) {
      p.diplo[l] <- p.haplo[i]*p.haplo[j]
      code.a[l] <- code.abc$a[i] + code.abc$a[j]
      code.b[l] <- code.abc$b[i] + code.abc$b[j]
      code.c[l] <- code.abc$c[i] + code.abc$c[j]
      l <- l + 1
    }
  
  N <- 1000
  N.diplo <- round(N*p.diplo)
  nind <- sum(N.diplo)
  
  geno.a <- rep(0, nind)
  geno.b <- rep(0, nind)
  geno.c <- rep(0, nind)
  
  
  l <- 1
  for (i in 1:n.diplo) {
    if (N.diplo[i]) {
      geno.a[l:(l+N.diplo[i]-1)] <- rep(code.a[i], N.diplo[i])
      geno.b[l:(l+N.diplo[i]-1)] <- rep(code.b[i], N.diplo[i])
      geno.c[l:(l+N.diplo[i]-1)] <- rep(code.c[i], N.diplo[i])
      l <- l + N.diplo[i]
    }
  }
  
  ## Genetic value
  gg.b <- eogCalc(pb, vb)
  gg.c <- eogCalc(pc, vc)
  
  ## Covariance between genetic means of B and C. var(pheno) = var(Gb) + var(Gc) + 2*cov(Gb, Gc)
  # Calculate the variance
  ggb <- gg.b[geno.b+1]
  ggc <- gg.c[geno.c+1]
  tovar <- vb + vc + 2*cov(ggb, ggc)
  if (!err)
    pheno <- sapply(1:nind, function(x) {gg.b[geno.b[x]+1] + gg.c[geno.c[x]+1]})
  else {
    pheno <- sapply(1:nind, function(x) {gg.b[geno.b[x]+1] + gg.c[geno.c[x]+1] + rnorm(1, 0, sqrt(tovar)*err)})# + rnorm(1)*sqrt(1-tovar)})
  }
  R2A <- summary(lm(pheno~geno.a))$r.squared
  
  #R2B <- summary(lm(pheno~geno.b))$r.squared
  #R2C <- summary(lm(pheno~geno.c))$r.squared
  #R2BC <- summary(lm(pheno~geno.b + geno.c))$r.squared
  #R2ABC <- summary(lm(pheno~geno.a + geno.b + geno.c))$r.squared
  ## Correlation between phenotype and genotype of A/C
  rBP <- cor(pheno, geno.b)
  rCP <- cor(pheno, geno.c)
  ## Correlation between genotype of B and genotype of A/C
  rAB <- cor(geno.b, geno.a)
  rBC <- cor(geno.b, geno.c)
  rAC <- cor(geno.a, geno.c)
  
  #Print(rBC)
  #res <- c(R2A, R2A.m)
  #  names(res) <- c("R2A", "R2B", "R2C", "R2AC", "R2ABC", "d2.AB", "d2.BC", "d2.AC",
  #                  "cor.A", "cor.C", "cor.a", "cor.c", "D.AB", "D.BC", "D.AC", "D.ABC")
  if(abs(abs(rBC) - 1) < 1e-3)
    R2A.m <- NA
  else 
    R2A.m <- calc_R2m_2(c(rAB, rAC), c(rBC, rBP, rCP))

  return(c(R2A, R2A.m, pa, pb, pc, rAB, rBC, rAC, rBP, rCP, 
           var(pheno), 2*cov(ggb, ggc)))

}

## Simulate three regulatory SNPs
simu3_sf <- function(vb, vc, vd, err = F) {
  
  hapcode.abc <- data.frame(a = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1), 
                            b = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1), 
                            c = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1), 
                            d = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1))
  
  #p.haplo <- c(pABCD,pABCd,paBCD,paBCd,pAbCD,pAbCd,pabCD,pabCd,pABcD,pABcd,paBcD,paBcd,pAbcD,pAbcd,pabcD,pabcd)
  p.haplo <- as.numeric(random_16())
  pa <- sum(p.haplo[as.logical(c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1))])
  pb <- sum(p.haplo[as.logical(c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1))])
  pc <- sum(p.haplo[as.logical(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1))])
  pd <- sum(p.haplo[as.logical(c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1))])
  
  n.haplo <- length(p.haplo)
  p.diplo <- rep(0, n.haplo^2)
  n.diplo <- length(p.diplo)
  dipcode.abc <- data.frame(a = rep(0, n.diplo), b = 0, c = 0, d = 0)
  
  l <- 1
  for (i in 1:n.haplo) 
    for (j in 1:n.haplo) {
      p.diplo[l] <- p.haplo[i]*p.haplo[j]
      dipcode.abc$a[l] <- hapcode.abc$a[i] + hapcode.abc$a[j]
      dipcode.abc$b[l] <- hapcode.abc$b[i] + hapcode.abc$b[j]
      dipcode.abc$c[l] <- hapcode.abc$c[i] + hapcode.abc$c[j]
      dipcode.abc$d[l] <- hapcode.abc$d[i] + hapcode.abc$d[j]
      l <- l + 1
    }
  
  N <- 1000
  N.diplo <- round(N*p.diplo)
  nind <- sum(N.diplo)
  geno.abc <- data.frame(a = rep(0, nind), b = 0, c = 0, d = 0)
  
  l <- 1
  for (i in 1:n.diplo) {
    if (N.diplo[i]) {
      geno.abc$a[l:(l+N.diplo[i]-1)] <- rep(dipcode.abc$a[i], N.diplo[i])
      geno.abc$b[l:(l+N.diplo[i]-1)] <- rep(dipcode.abc$b[i], N.diplo[i])
      geno.abc$c[l:(l+N.diplo[i]-1)] <- rep(dipcode.abc$c[i], N.diplo[i])
      geno.abc$d[l:(l+N.diplo[i]-1)] <- rep(dipcode.abc$d[i], N.diplo[i])
      l <- l + N.diplo[i]
    }
  }
  
  gg.b <- eogCalc(pb, vb)
  gg.c <- eogCalc(pc, vc)
  gg.d <- eogCalc(pd, vd)
 
  ggb <- gg.b[geno.abc$b+1]
  ggc <- gg.c[geno.abc$c+1]
  ggd <- gg.d[geno.abc$d+1]
  tovar <- vb + vc + vd + 2*(cov(ggb, ggc) + cov(ggc, ggd) + cov(ggb, ggd))
  if (!err)
    pheno <- sapply(1:nind, function(x) {gg.b[geno.abc$b[x]+1] + gg.c[geno.abc$c[x]+1] + gg.d[geno.abc$d[x]+1]})
  else {
    pheno <- sapply(1:nind, function(x) {gg.b[geno.abc$b[x]+1] + gg.c[geno.abc$c[x]+1] + gg.d[geno.abc$d[x]+1] + rnorm(1, 0, sqrt(tovar)*err)})
  }
  R2A <- summary(lm(pheno~geno.abc$a))$r.squared
  
  ## Correlation between phenotype and genotype of A/C
  rBP <- cor(pheno, geno.abc$b)
  rCP <- cor(pheno, geno.abc$c)
  rDP <- cor(pheno, geno.abc$d)
  
  ## Correlation between genotype of B and genotype of A/C
  rAB <- cor(geno.abc$a, geno.abc$b)
  rBC <- cor(geno.abc$b, geno.abc$c)
  rAC <- cor(geno.abc$a, geno.abc$c)
  rAD <- cor(geno.abc$a, geno.abc$d)
  rBD <- cor(geno.abc$b, geno.abc$d)
  rCD <- cor(geno.abc$c, geno.abc$d)

  # SNP B, C, D should be not linked
  if(abs(abs(rBC) - 1) < 1e-3 | abs(abs(rBD) - 1) < 1e-3 | abs(abs(rCD) - 1) < 1e-3)
    R2A.m <- NA
  else 
    R2A.m <- calc_R2m_3(c(rAB, rAC, rAD), c(rBC, rBD, rCD, rBP, rCP, rDP))
  return(c(R2A, R2A.m, pa, pb, pc, pd, rAB, rBC, rAC, rAD, rBD, rCD, rBP, rCP, rDP))
  
}

## Simulate four regulatory SNPs
simu4_sf <- function(vb, vc, vd, ve, err = F) {
  
  hapcode.abc <- data.frame(a = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1), 
                            b = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1), 
                            c = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1), 
                            d = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
                            e = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
  
  #p.haplo <- c(pABCD,pABCd,paBCD,paBCd,pAbCD,pAbCd,pabCD,pabCd,pABcD,pABcd,paBcD,paBcd,pAbcD,pAbcd,pabcD,pabcd)
  p.haplo <- as.numeric(random_32())
  pa <- sum(p.haplo[as.logical(c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1))])
  pb <- sum(p.haplo[as.logical(c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1))])
  pc <- sum(p.haplo[as.logical(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1))])
  pd <- sum(p.haplo[as.logical(c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1))])
  pe <- sum(p.haplo[as.logical(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))])
  
  n.haplo <- length(p.haplo)
  p.diplo <- rep(0, n.haplo^2)
  n.diplo <- length(p.diplo)
  dipcode.abc <- data.frame(a = rep(0, n.diplo), b = 0, c = 0, d = 0, e=0)
  
  l <- 1
  for (i in 1:n.haplo) 
    for (j in 1:n.haplo) {
      p.diplo[l] <- p.haplo[i]*p.haplo[j]
      dipcode.abc$a[l] <- hapcode.abc$a[i] + hapcode.abc$a[j]
      dipcode.abc$b[l] <- hapcode.abc$b[i] + hapcode.abc$b[j]
      dipcode.abc$c[l] <- hapcode.abc$c[i] + hapcode.abc$c[j]
      dipcode.abc$d[l] <- hapcode.abc$d[i] + hapcode.abc$d[j]
      dipcode.abc$e[l] <- hapcode.abc$e[i] + hapcode.abc$e[j]
      l <- l + 1
    }
  
  N <- 1000
  N.diplo <- round(N*p.diplo)
  nind <- sum(N.diplo)
  geno.abc <- data.frame(a = rep(0, nind), b = 0, c = 0, d = 0, e=0)
  
  l <- 1
  for (i in 1:n.diplo) {
    if (N.diplo[i]) {
      geno.abc$a[l:(l+N.diplo[i]-1)] <- rep(dipcode.abc$a[i], N.diplo[i])
      geno.abc$b[l:(l+N.diplo[i]-1)] <- rep(dipcode.abc$b[i], N.diplo[i])
      geno.abc$c[l:(l+N.diplo[i]-1)] <- rep(dipcode.abc$c[i], N.diplo[i])
      geno.abc$d[l:(l+N.diplo[i]-1)] <- rep(dipcode.abc$d[i], N.diplo[i])
      geno.abc$e[l:(l+N.diplo[i]-1)] <- rep(dipcode.abc$e[i], N.diplo[i])
      l <- l + N.diplo[i]
    }
  }
  
  gg.b <- eogCalc(pb, vb)
  gg.c <- eogCalc(pc, vc)
  gg.d <- eogCalc(pd, vd)
  gg.e <- eogCalc(pe, ve)
  
  ggb <- gg.b[geno.abc$b+1]
  ggc <- gg.c[geno.abc$c+1]
  ggd <- gg.d[geno.abc$d+1]
  gge <- gg.e[geno.abc$e+1]
  
  tovar <- vb + vc + vd + ve + 
           2*(cov(ggb, ggc) + cov(ggc, ggd) + cov(ggb, ggd)) + 
           2*(cov(ggb, gge) + cov(ggc, gge) + cov(ggd, gge))
  if (!err)
    pheno <- sapply(1:nind, function(x) {gg.b[geno.abc$b[x]+1] + gg.c[geno.abc$c[x]+1] + gg.d[geno.abc$d[x]+1] + gg.e[geno.abc$e[x]+1]})
  else {
    pheno <- sapply(1:nind, function(x) {gg.b[geno.abc$b[x]+1] + gg.c[geno.abc$c[x]+1] + gg.d[geno.abc$d[x]+1] + gg.e[geno.abc$e[x]+1] + rnorm(1, 0, sqrt(tovar)*err)})
  }
  R2A <- summary(lm(pheno~geno.abc$a))$r.squared
  
  ## Correlation between phenotype and genotype of A/C
  rBP <- cor(pheno, geno.abc$b)
  rCP <- cor(pheno, geno.abc$c)
  rDP <- cor(pheno, geno.abc$d)
  rEP <- cor(pheno, geno.abc$e)
 
  ## Correlation between genotype of B and genotype of A/C
  rAB <- cor(geno.abc$a, geno.abc$b)
  rBC <- cor(geno.abc$b, geno.abc$c)
  rAC <- cor(geno.abc$a, geno.abc$c)
  rAD <- cor(geno.abc$a, geno.abc$d)
  rBD <- cor(geno.abc$b, geno.abc$d)
  rCD <- cor(geno.abc$c, geno.abc$d)
  rAE <- cor(geno.abc$a, geno.abc$e)
  rBE <- cor(geno.abc$b, geno.abc$e)
  rCE <- cor(geno.abc$c, geno.abc$e)
  rDE <- cor(geno.abc$d, geno.abc$e)

  # SNP B, C, D, E should be not linked
  if(abs(abs(rBC) - 1) < 1e-3 | abs(abs(rBD) - 1) < 1e-3 | abs(abs(rCD) - 1) < 1e-3 | abs(abs(rBE) - 1) < 1e-3 | abs(abs(rCE) - 1) < 1e-3 | abs(abs(rDE) - 1) < 1e-3)
    R2A.m <- NA
  else 
    R2A.m <- calc_R2m_4(c(rAB, rAC, rAD, rAE), c(rBC, rBD, rBE, rCD, rCE, rDE, rBP, rCP, rDP, rEP))
  return(c(R2A, R2A.m))#, pa, pb, pc, pd, rAB, rBC, rAC, rAD, rBD, rCD, rBP, rCP, rDP))
  
}

## Calculate genetic mean using the method from PLINK by Shaun Purcell
eogCalc <- function (pa, v, dom = 0) {
  if (pa > 0.5)
    pa <- 1 - pa
  if (pa < 1e-8)
    stop("MAF should be > 0")
  p <- pa
  q <-  1-p
  # to assign the A first, them calculate the va
  A <-  sqrt( ( v ) / ( (2*p*q)* (1+dom*(q-p))*(1+dom*(q-p)) + (2*p*q*dom)*(2*p*q*dom) ) )
  D <-  dom * A
  gBB <-  A  - (A*(p-q) + (2*p*q*D))
  gAB <-  D  - (A*(p-q) + (2*p*q*D))
  gAA <- -A  - (A*(p-q) + (2*p*q*D))
  return(c(gAA, gAB, gBB)) 
}

## The random_x functions generate x random numbers that sum to 1
## This is used to mimic the haplotype frequencies that sum to 1
random_4 <- function() {
  pa <- pb <- 0
  
  ## make sure the MAF is within legal range
  while (pa < 0.05 | pb < 0.05 | pa > 0.5 | pb > 0.5) {
    n.zero <- sample(0:2, 1) ## introduce 0 values to mimic real haptypes that ususally have multiple zeros
    tmp <- sample(1:100, 4-n.zero) ## otherwise won't give R2 to 1, this happens when some haplotypes are zeros 
    tmp <- c(tmp, rep(0, n.zero))
    tmp <- sample(tmp/sum(tmp))
    pa <- sum(tmp[as.logical(c(0, 0, 1, 1))])
    pb <- sum(tmp[as.logical(c(0, 1, 0, 1))])
  }
  return(tmp)
}

random_8 <- function() {
  pa <- pb <- pc <- 0
 
  ## make sure the MAF is within legal range  | pa > 0.5 | pb > 0.5 | pc > 0.5
  while (pa < 0.05 | pb < 0.05 | pc < 0.05| pa > 0.95 | pb > 0.95 | pc > 0.95) {
    n.zero <- sample(0:6, 1) ## introduce 0 values to mimic real haptypes that ususally have multiple zeros
    tmp <- sample(1:100, 8-n.zero, replace = T)###-n.zero otherwise won't give R2 to 1, this happens when some haplotypes are zeros 
    tmp <- c(tmp, rep(0, n.zero))
    tmp <- tmp/sum(tmp)
    pa <- sum(tmp[as.logical(c(0, 0, 1, 1, 0, 0, 1, 1))])
    pb <- sum(tmp[as.logical(c(0, 1, 0, 1, 0, 1, 0, 1))])
    pc <- sum(tmp[as.logical(c(0, 1, 1, 0, 1, 0, 0, 1))])
  }
  return(tmp)
}

random_16 <- function() {
  pa <- pb <- pc <- pd <- 0
  
  ## Make sure the MAF is within legal range
  while (pa < 0.05 | pb < 0.05 | pc < 0.05 | pd < 0.05 | pa > 0.5 | pb > 0.5 | pc > 0.5 | pd > 0.5) {
    n.zero <- sample(0:14, 1) ## introduce 0 values to mimic real haptypes that ususally have multiple zeros
    tmp <- sample(1:100, 16-n.zero) ## otherwise won't give R2 to 1, this happens when some haplotypes are zeros 
    tmp <- c(tmp, rep(0, n.zero))
    tmp <- sample(tmp/sum(tmp))
    pa <- sum(tmp[as.logical(c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1))])
    pb <- sum(tmp[as.logical(c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1))])
    pc <- sum(tmp[as.logical(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1))])
    pd <- sum(tmp[as.logical(c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1))])
  }
  return(tmp)
}

random_32 <- function() {
  pa <- pb <- pc <- pd <- pe <- 0
  
  ## Make sure the MAF is within legal range
  while (pa < 0.05 | pb < 0.05 | pc < 0.05 | pd < 0.05 | pe <0.05 | pa > 0.5 | pb > 0.5 | pc > 0.5 | pd > 0.5 | pe > 0.5) {
    n.zero <- sample(0:30, 1) ## introduce 0 values to mimic real haptypes that ususally have multiple zeros
    tmp <- sample(1:100, 32-n.zero) ## otherwise won't give R2 to 1, this happens when some haplotypes are zeros 
    tmp <- c(tmp, rep(0, n.zero))
    tmp <- sample(tmp/sum(tmp))
    pa <- sum(tmp[as.logical(c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1))])
    pb <- sum(tmp[as.logical(c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1))])
    pc <- sum(tmp[as.logical(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1))])
    pd <- sum(tmp[as.logical(c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1))])
    pe <- sum(tmp[as.logical(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))])
  }
  return(tmp)
}
