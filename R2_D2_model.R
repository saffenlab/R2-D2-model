require(MASS)
require(RcppEigen)
require(data.table)

getPhenotypeBC <- function (gene, probe.name, expr.dir, probe.dir) {
  expr <- fread(expr.dir, head = T)
  probe <- fread(probe.dir, head = T)
  setkey(probe, "Gene_Symbol")
  probeName <- probe[gene]$ID
  if (is.na(match(probe.name, probeName)))
    stop(paste0("Cannot find probe ", probe.name, " for gene ", gene))
  else
    expression <- as.numeric(unlist(expr[, probe.name, with = F]))
  write.table(cbind(expr[, c("FID", "IID"), with = F], expression), paste0(gene, "_", probe.name, "_phenotype-BC.txt"), row.names = F, quote = F, sep = "\t")
}

getPhenotypeFCTX <- function (gene, probe.name, expr.dir, probe.dir) {
  expr <- fread(expr.dir, head = T)
  probe <- fread(probe.dir, head = T)
  setkey(probe, "Gene_Symbol")
  probeName <- probe[gene]$ID
  if (is.na(match(probe.name, probeName)))
    stop(paste0("Cannot find probe ", probe.name, " for gene ", gene))
  else
    expression <- as.numeric(unlist(expr[, probe.name, with = F]))
  
	write.table(cbind(expr[, c("FID", "IID"), with = F], expression), paste0(gene, "_", probe.name, "_phenotype-FCTX.txt"), row.names = F, quote = F, sep = "\t")
}

getPhenotypeTCTX <- function (gene, probe.name, expr.dir, probe.dir) {
  expr <- fread(expr.dir, head = T)
  probe <- fread(probe.dir, head = T)
  setkey(probe, "Gene_Symbol")
  probeName <- probe[gene]$ID
  if (is.na(match(probe.name, probeName)))
    stop(paste0("Cannot find probe ", probe.name, " for gene ", gene))
  else
    expression <- as.numeric(unlist(expr[, probe.name, with = F]))
  
  write.table(cbind(expr[, c("FID", "IID"), with = F], expression), paste0(gene, "_", probe.name, "_phenotype-TCTX.txt"), row.names = F, quote = F, sep = "\t")
}

getPhenotypeCERE <- function (gene, probe.name, expr.dir, probe.dir) {
  expr <- fread(expr.dir, head = T)
  probe <- fread(probe.dir, head = T)
  setkey(probe, "Gene_Symbol")
  probeName <- probe[gene]$ID
  if (is.na(match(probe.name, probeName)))
    stop(paste0("Cannot find probe ", probe.name, " for gene ", gene))
  else
    expression <- as.numeric(unlist(expr[, probe.name, with = F]))
  
  write.table(cbind(expr[, c("FID", "IID"), with = F], expression), paste0(gene, "_", probe.name, "_phenotype-CERE.txt"), row.names = F, quote = F, sep = "\t")
}

getPhenotypePONS <- function (gene, probe.name, expr.dir, probe.dir) {
  expr <- fread(expr.dir, head = T)
  probe <- fread(probe.dir, head = T)
  setkey(probe, "Gene_Symbol")
  probeName <- probe[gene]$ID
  if (is.na(match(probe.name, probeName)))
    stop(paste0("Cannot find probe ", probe.name, " for gene ", gene))
  else
    expression <- as.numeric(unlist(expr[, probe.name, with = F]))
  
  write.table(cbind(expr[, c("FID", "IID"), with = F], expression), paste0(gene, "_", probe.name, "_phenotype-PONS.txt"), row.names = F, quote = F, sep = "\t")
}

pBar <- function (x, percent = 0) {
  cat("\f")
  if (x + percent < 10)
    print(noquote(paste0("  ", x + percent, "%")))
  else
    print(noquote(paste0(" ", x + percent, "%")))
}

## Main function
comTest_matrix <- function (geno, qassoc, num, save_interm=FALSE) {
  if (num < 0) {stop("num must be non-negative!")}
  start <- Sys.time()
  pBar(0)
  ## input data
  R2B.regss <- extract_R2(qassoc)
  P <- extract_P(qassoc)
  if (num <= 4)
    p.05 <- (P < 0.05) ## choose only significant SNPs
  else
    p.05 <- (P < 0.0001) ## 
  pBar(1)
  rGP <- calc_rGP(qassoc)
  pBar(2)
  rGG <- calc_rGG(geno)
  rGG.05 <- as.matrix(rGG[, p.05]) ## choose only significant SNPs
  pBar(3)
  
  snpName <- colnames(rGG)
  rGG_uniq <- make_uniq(rGG.05)
  snpName_uniq <- colnames(rGG_uniq)
  nSNPs <- length(snpName_uniq)
  pBar(4)
  
  if (num > 0) {
  ## Combination calculation
  N <- choose(nSNPs, num)
  combi <- combn(nSNPs, num)
  res <- data.frame(SNP = "a", NRMSE = 0, R2.adj = rep(0.0, N), AIC = 0, BIC = 0, R2.T = 0,
                    R2B = 0, R2C = 0, R2D = 0, R2E = 0, stringsAsFactors = F)
  cnt <- 1

  NRMSE.min <- Inf
  n <- length(R2B.regss)
  
  pBar(5)
  for (m in 1:N) {
    if (floor(m*93/N) - floor((m-1)*93/N) == 1)
      pBar(5, floor(m*93/N))
    
    ori.index <- match(snpName_uniq[combi[, m]], snpName)
    if (num == 1) {
      R2m <- R2B.regss[ori.index]*(rGG[, ori.index])^2
      R2t <- R2B.regss[ori.index]
    } else if (num == 2){
      i <- ori.index[1]
      j <- ori.index[2]
      tmp <- c(rGG[i, j], rGP[i], rGP[j])
      R2m <- apply(rGG[, c(i, j)], 1, function(x){calc_R2m_2(x, tmp)})
      R2t <- calc_R2t(rGG[c(i, j), c(i, j)], rGP[c(i, j)])
    } else if (num == 3) {
      i <- ori.index[1]
      j <- ori.index[2]
      k <- ori.index[3]
      tmp.six <- c(rGG[i, j], rGG[i, k], rGG[j, k], rGP[i], rGP[j], rGP[k])
      R2m <- apply(rGG[, c(i, j, k)], 1, function(x){calc_R2m_3(x, tmp.six)})
      R2t <- calc_R2t(rGG[c(i, j, k), c(i, j, k)], rGP[c(i, j, k)])
    } else if (num == 4) {
      i <- ori.index[1]
      j <- ori.index[2]
      k <- ori.index[3]
      l <- ori.index[4]
      tmp.six <- c(rGG[i, j], rGG[i, k], rGG[i, l], rGG[j, k], rGG[j, l], rGG[k, l], rGP[i], rGP[j], rGP[k], rGP[l])
      R2m <- apply(rGG[, c(i, j, k, l)], 1, function(x){calc_R2m_4(x, tmp.six)})
      R2t <- calc_R2t(rGG[c(i, j, k, l), c(i, j, k, l)], rGP[c(i, j, k, l)])
    }
    else {stop("Can handle combinations of 4 SNPs at most")}
    
    lm_stats <- fastStats(R2m, R2B.regss)
    R2.adj <- lm_stats[1]
    R2.adj <- ifelse(is.na(R2.adj), 0, R2.adj)
    
    ## Minimize errors for the major peaks only to avoid selecting SNPs with "small" contributions to the variance of mRNA expression.
    index.major <- (R2B.regss > max(R2B.regss)/3)
    R2m.major <- R2m[index.major]
    R2B.regss.major <- R2B.regss[index.major]
    NRMSE <- sqrt(mean((R2m.major-R2B.regss.major)^2))/(max(R2B.regss.major)-min(R2B.regss.major))*100
    ## continue the program when NRMSE is infinite
    if (NRMSE == 0 || is.infinite(NRMSE) || is.na(NRMSE) || is.na(R2t) || is.na(R2.adj) || (sum(is.na(R2m)) > length(R2m)/3) || (sum(is.infinite(R2m)) > length(R2m)/3)) {
      next
    }

    # Output AIC/BIC to file
    aic <- lm_stats[2]
    bic <- lm_stats[3]
    
    if (num == 4) {
      res[cnt, ] <- list(paste(snpName[ori.index], collapse = "\t"), round(NRMSE, digits = 6),
                       round(R2.adj, digits = 6), round(aic, digits = 6), round(bic, digits = 6),
                       round(R2t, digits = 6), R2B.regss[i], R2B.regss[j], R2B.regss[k], R2B.regss[l])
    } else if (num == 3) {
      res[cnt, ] <- list(paste(snpName[ori.index], collapse = "\t"), round(NRMSE, digits = 6),
                       round(R2.adj, digits = 6), round(aic, digits = 6), round(bic, digits = 6),
                       round(R2t, digits = 6), R2B.regss[i], R2B.regss[j], R2B.regss[k], NA)
    } else if (num == 2) {
      res[cnt, ] <- list(paste(snpName[ori.index], collapse = "\t"), round(NRMSE, digits = 6),
                       round(R2.adj, digits = 6), round(aic, digits = 6), round(bic, digits = 6),
                       round(R2t, digits = 6), R2B.regss[i], R2B.regss[j], NA, NA)
    } else {
      res[cnt, ] <- list(paste(snpName[ori.index], collapse = "\t"), round(NRMSE, digits = 6),
                       round(R2.adj, digits = 6), round(aic, digits = 6), round(bic, digits = 6),
                       round(R2t, digits = 6), R2B.regss[ori.index], NA, NA, NA)
    }
    
    if (NRMSE.min > NRMSE) {
      NRMSE.min <- NRMSE 
      R2_predicted <- R2m
      snp.best <- ori.index
    }
    cnt <- cnt + 1
  }
  res <- res[1:(cnt - 1), ]
  res <- res[order(res$NRMSE, decreasing = F) , ]
  
  pBar(98)
  rGG.best.snps <- as.data.frame(rGG[, snp.best])
  colnames(rGG.best.snps) <- snpName[snp.best]
  plot.file <- cbind(data.frame(SNP = snpName, R2 = R2B.regss), rGG.best.snps, P, R2_predicted)
  
  # Output file
  pBar(99)
  colnames(res)[1] <- paste0(paste0('SNP_', seq(1, num)), collapse = "\t")
  write.table(res, paste0(qassoc, "-", num, "snps-m.txt"), quote=F, sep="\t", row.names=F) 
  write.table(plot.file, paste0(qassoc, "-", num, "snps-plot-m.txt"), quote=F, sep="\t", row.names=F) 
  
  } ## <-- end of if (num > 0)
  else {
    write.table(data.frame(SNP = snpName, R2 = R2B.regss, P=P), paste0(qassoc, "-", num, "snps-plot-m.txt"), quote=F, sep="\t", row.names=F) 
  }
  ## Save intermediate files
  if (save_interm) {
    snp_name_df = cbind(snpName, rGG)
    colnames(snp_name_df) <- c('rsID', snpName)
    write.table(snp_name_df, paste0(qassoc, "-", num, "snps-r2GG.txt"), quote=F, sep="\t", row.names=F) 
  }
  
  end <- Sys.time() 
  pBar(100)
  print(noquote(paste0("Run time: ", round(as.numeric(end - start, units='secs'), digits = 2), "s")))
}

## Calculate R2m  for four regulatory variants (R2m = "predicted" R2; m = matrix.)
calc_R2m_4 <- function(r.Ax, tmp.six) {
  
  a12 <- r.Ax[1]
  a13 <- r.Ax[2]
  a14 <- r.Ax[3]
  a15 <- r.Ax[4]
  
  a23 <- tmp.six[1]
  a24 <- tmp.six[2]
  a25 <- tmp.six[3]
  
  a34 <- tmp.six[4]
  a35 <- tmp.six[5]
  
  a45 <- tmp.six[6]
  
  B <- tmp.six[7]
  C <- tmp.six[8]
  D <- tmp.six[9]
  E <- tmp.six[10]
  
  a <- matrix(c(1, a12, a13, a14, a15, 
                a12, 1, a23, a24, a25, 
                a13, a23, 1, a34, a35, 
                a14, a24, a34, 1, a45, 
                a15, a25, a35, a45, 1), nrow = 5, ncol=5, byrow = T)
  
  # Calculate the adjugate matrix (identified as "adjunctMatrix" within the script below)
  minor <- function(A, i, j) det(A[-i, -j])
  cofactor <- function(A, i, j) (-1)^(i+j) * minor(A, i, j)
  
  adjunctMatrix <- function(A) {
    n <- nrow(A)
    t(outer(1:n, 1:n, Vectorize(function(i, j) cofactor(A, i, j)
  )))}
  
  b = adjunctMatrix(a)
  
  if (any(is.na(as.numeric(b)))) {
    R2m = 0
  } else {
    if (abs(b[1, 1]) > 1e-7)
      R2m <- (sum(c(B, C, D, E)*c(b[1, 2], b[1, 3], b[1, 4], b[1, 5]))/b[1, 1])^2
    else
      R2m <- Inf
  }
  return(R2m)
}

## Calculate R2m for three regualtory variants
calc_R2m_3 <- function(r.Ax, tmp.six) {
  
  a <- r.Ax[1]
  c <- r.Ax[2]
  f <- r.Ax[3]
  b <- tmp.six[1]
  d <- tmp.six[2]
  e <- tmp.six[3]
  
  b11 <- 1 + 2*b*d*e - b^2 - d^2 - e^2
  b12 <- b*c + d*f + a*e^2 - a - c*d*e - f*b*e
  b13 <- a*b + c*d^2 + e*f - a*d*e - c - b*d*f
  b14 <- a*d + c*e + b^2*f - a*b*e - f - b*c*d
  
  if (abs(b11) > 1e-7)
    R2m <- (sum(tmp.six[4:6]*c(b12, b13, b14))/b11)^2
  else
    R2m <- Inf
  return(R2m)
}

## Calculate R2m for two regulatory variants
calc_R2m_2 <- function(r.Ax, tmp) {
  
  a <- r.Ax[1]
  c <- r.Ax[2]
  b <- tmp[1]
  rBP <- tmp[2]
  rCP <- tmp[3]
  
  if (abs(abs(b) - 1) > 1e-7)
    R2m <- (1/(1-b^2)^2)*(rBP^2*(a-c*b)^2 + rCP^2*(c-a*b)^2 + 2*(a-c*b)*(c-a*b)*rBP*rCP)
  else
    R2m <- Inf
  return(R2m)
}

extract_R2 <- function (fileName) {
  qassoc <- read.table(paste0(fileName, ".qassoc"), header = T)
  return(qassoc$R2)
}

extract_P <- function (fileName) {
  qassoc <- read.table(paste0(fileName, ".qassoc"), header = T)
  return(qassoc$P)
}

calc_rGP <- function (fileName) {
  qassoc <- read.table(paste0(fileName, ".qassoc"), header = T)
  return(sign(qassoc$BETA)*sqrt(qassoc$R2))
}

calc_rGG <- function (fileName) {
  geno <- read.table(paste0(fileName, ".raw"), header = T)
  geno <- geno[order(geno$FID), 7:ncol(geno)]
  nSNPs <- ncol(geno)
  
  snpName <- as.character(sapply(colnames(geno), function(x){unlist(strsplit(x, "_"))[1]}))
  r.matrix <- matrix(0, nr = nSNPs, nc = nSNPs, dimnames = list(snpName, snpName))
  
  for (i in 1:nSNPs)
    for (j in 1:nSNPs) {
      if (i <= j)
        r.matrix[j, i] <- round(cor(geno[, i], geno[, j], use = "complete.obs"), digits = 4)
      else
        r.matrix[j, i] <-  r.matrix[i, j]
    }
  return(r.matrix)
}

make_uniq <- function (matrix) {
  return(matrix[, !duplicated(t(matrix))])
}

# Calculate R2 for specified SNPs 
calc_R2 <- function(geno, qassoc, snplist) {
  R2B.regss <- extract_R2(qassoc)
  rGG <- calc_rGG(geno)
  snpName <- colnames(rGG)
  rGG.best.snps <- as.data.frame(rGG[, snplist]^2)
  colnames(rGG.best.snps) <- snplist
  df <- cbind(R2B.regss, rGG.best.snps)
  lm1 <- lm(formula(paste0("R2B.regss~",paste(snplist,collapse="+"))), data = df)
  sumLm <- summary(lm1)
  adjusted_R2 <- round(sumLm$adj.r.squared, digits = 8)
  P <- format(1-pf(sumLm$fstat[1], sumLm$fstat[2], sumLm$fstat[3]), scientific = T, digits=10)
  return(c(adjusted_R2, P))
}

calc_R2t <- function(rgg, rgp) {
  R2t = tryCatch({abs(t(rgp)%*%solve(rgg)%*%rgp)}, error=function(e) {return(NA)}, finally = {})
  return(R2t)
}

fastStats <- function(X, Y) {
  
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  X.n <- nrow(X)
  X.p <- ncol(X)
  res <- fastLm(cbind(1, X), Y, method = 2)
  mean_Y <- mean(Y)
  R2 <- sum((res$fit-mean_Y)^2)/sum((Y-mean_Y)^2)
  R2.adjusted <- 1 - (X.n-1)/(X.n-X.p-1)*(1-R2)
  
  # Formulas are from Wikipedia
  rss = sum((res$fit-Y)^2)
  aic = 2*X.p + X.n*log(rss)
  bic = X.p*log(X.n) + X.n*log(rss/X.n)
  
  return(c(R2.adjusted, aic, bic))
}

familyPlot_matrix <- function (plot_file, title, color = c('black', 'red', 'green', 'brown')) {
  
  E.all <- read.table(paste0(plot_file, ".txt"), header=T, quote="\\")
  
  plotTitle <- paste0(title, "_Matrix_Model")
  
  fileName <- paste0(plotTitle, ".pdf") 
  pdf(fileName, height=100, width=300) # 90-500
  
  nSNPs <- nrow(E.all) # total number of snps
  snpName <- E.all$SNP
  
  fill_colors <- rep("red", nSNPs)
  
  for(i in 1:nSNPs)
    fill_colors[i] <- ifelse(E.all$P[i] < 0.05, "#99CCFF", "#cccccc")
  
  ## if no iSNP provided, only plot the R2 as barplot
  if (ncol(E.all) == 3) {
    par(mfrow = c(1, 1), mar=c(145,30,145,2), las = 2, cex.lab=15, cex.axis=15)#, mgp=c(40, 1, 0))#
    mp <- barplot(E.all$R2, names.arg=snpName,
                  col=fill_colors, horiz=FALSE, border=NA, xpd = F,cex.names=8)
  } 
  else 
  {
  par(mfrow = c(2, 1), mar=c(40,30,5,2), las = 2, cex.lab=15, cex.axis=15)#, mgp=c(40, 1, 0))#
  
  snps <- colnames(E.all)[-c(1, 2, ncol(E.all)-1, ncol(E.all))]#
  m <- length(snps) # number of index snps
  
  lm_stats <- round(fastStats(E.all$R2, E.all$R2_predicted), digits = 3)
  R2.adj <- lm_stats[1]
  #index.major <- (E.all$R2 >= min(E.all$R2) )# max(E.all$R2)/3)
  #R2m.major <- E.all$R2_predicted[index.major]
  #R2B.regss.major <- E.all$R2[index.major]
  #NRMSE <- sqrt(mean((R2m.major-R2B.regss.major)^2))/(max(R2B.regss.major)-min(R2B.regss.major))*100
  #NRMSE <- round(NRMSE, digits = 2)
  #NRMSE <- round(sqrt(mean((E.all$R2_predicted-E.all$R2)^2))/(max(E.all$R2)-min(E.all$R2))*100, digits = 2)
  #Plot predicted R2 vs. estimated R2 ("estimated" R2 = R2 values obtained from linear regression analysis of mRNA expression vs SNP genotype,
  #calculated independently for each SNP in the chromosome ROI dataset).
  mp <- barplot(E.all$R2, names.arg=snpName,
                col=fill_colors, horiz=FALSE, border=NA, xpd = F,cex.names=8)
  lines(mp, E.all$R2_predicted, col="darkblue", lwd=50) #blue2
  legend("topright",legend=paste0("predicted ", paste0(snps, collapse = " - ")),
         lty=1,lwd=50,col="darkblue", bty="n", cex=15, title=plotTitle) 
  legend("topleft",legend=as.expression(bquote("adjusted " ~ R[Model]^2 ~ " = " ~ .(R2.adj))), bty="n", cex=15) 
  
  #Plot index snp lines
  #Color <- c('red', 'black', 'green')
  #Color <- c('blue', 'purple', 'yellow') 
  mp <- barplot(E.all$R2, names.arg=snpName,
                col=fill_colors, horiz=FALSE, border=NA, xpd = F,cex.names=8)
  for (i in 1:m) {
    lines(mp, E.all[, i+2]^2*E.all$R2[which(snps[i] == snpName)], col = color[i], lwd = 50)
  }
  
  
  legend("topright",legend=snps,
         lty=1,lwd=50,col=color, bty="n",cex=15,title=plotTitle) 
  #Label for index SNPs on the bar
  for (i in 1:m) {
    label_all <- rep(NA, nSNPs)
    index_i <- which(snps[i] == snpName)
    label_all[index_i] <- i#c(3, 1, 2)[i]
    text(mp, E.all$R2[index_i]/2, labels = label_all, cex = 20, col = color[i])
  }
  
  }  ## end of else
  
  dev.off()
}

familyPlot_matrix_snps <- function (geno, qassoc, snplist, color=c('black', 'red', 'green', 'brown')) {
  
  R2B.regss <- extract_R2(qassoc)
  P <- extract_P(qassoc)
  rGP <- calc_rGP(qassoc)
  rGG <- calc_rGG(geno)
  
  snpName <- colnames(rGG)
  ori.index <- match(snplist, snpName)
  num <- length(snplist)
  if (num == 1) {
    R2m <- R2B.regss[ori.index]*(rGG[, ori.index])^2
  } else if (num == 2){
    i <- ori.index[1]
    j <- ori.index[2]
    tmp <- c(rGG[i, j], rGP[i], rGP[j])
    R2m <- apply(rGG[, c(i, j)], 1, function(x){calc_R2m_2(x, tmp)})
  } else if (num == 3) {
    i <- ori.index[1]
    j <- ori.index[2]
    k <- ori.index[3]
    tmp.six <- c(rGG[i, j], rGG[i, k], rGG[j, k], rGP[i], rGP[j], rGP[k])
    R2m <- apply(rGG[, c(i, j, k)], 1, function(x){calc_R2m_3(x, tmp.six)})
  } else if (num == 4) {
    i <- ori.index[1]
    j <- ori.index[2]
    k <- ori.index[3]
    l <- ori.index[4]
    tmp.six <- c(rGG[i, j], rGG[i, k], rGG[i, l], rGG[j, k], rGG[j, l], rGG[k, l], rGP[i], rGP[j], rGP[k], rGP[l])
    R2m <- apply(rGG[, c(i, j, k, l)], 1, function(x){calc_R2m_4(x, tmp.six)})
  }
  else {stop("Can handle combinations of 3 SNPs at most")}
  
  rGG.best.snps <- as.data.frame(rGG[, snplist])
  colnames(rGG.best.snps) <- snplist
  
  R2_predicted <- R2m
  plot.data <- cbind(data.frame(SNP = snpName, R2 = R2B.regss), rGG.best.snps, P, R2_predicted)
  plot.file <- paste0(qassoc, "-", paste0(snplist, collapse = "_"), "-plot-m")
  write.table(plot.data, paste0(plot.file, ".txt"), quote=F, sep="\t", row.names=F) 
  
  title=paste0(c(geno, snplist), collapse = "_")
  familyPlot_matrix(plot.file, title, color=color)
  if (file.exists(paste0(plot.file, ".txt")))
    file.remove(paste0(plot.file, ".txt"))
}

findSNPFamily <- function(plot.file, out.file) {
  r2 <- read.table(plot.file, header = T, stringsAsFactors = F)
  index.snps <- 3:(which(colnames(r2) == 'P') - 1)
  r2.index.snps <- r2[, index.snps]
  family.index <- apply(r2.index.snps, 1, function(x){x[which.max(x^2)] <- 1; not.max <- (index.snps-2)[-which.max(x^2)]; x[not.max] <- 0; return(x)})
  r2[, index.snps] <- t(family.index)
  write.table(r2[, c(1, 2, index.snps)], out.file, quote = F, sep = '\t', row.names = F)
}