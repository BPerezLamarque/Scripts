
#### This script contains amended functions of ParaFit and PACo from the R-packages ape (Paradis et al., 2004) and paco (Hutchinson et al., 2017))
#### The following functions have been amended to handle the cases were the number of microbial strains is low (<4), which can generate bugs in the original functions during the permutations. 


#### Contact: Benoît Perez-Lamarque (benoit.perez.lamarque@gmail.com)

#### For more details about this script, see "Comparing different computational approaches for detecting long-term vertical transmission in host-associated microbiota", Perez-Lamarque B & Morlon H, in prep.


parafit_test <- function (host.D, para.D, HP, nperm = 999, test.links = FALSE, 
                          seed = NULL, correction = "none", silent = FALSE,
                          nullmodel=1) {
  epsilon <- sqrt(.Machine$double.eps)
  if (is.null(seed)) {
    runif(1)
    seed <- .Random.seed[trunc(runif(1, 1, 626))]
  }
  HP <- as.matrix(HP)
  host.D <- as.matrix(host.D)
  host.pc <- pcoa_test(host.D, correction = correction)
  if (host.pc$correction[2] == 1) {
    if (min(host.pc$values[, 2]) < -epsilon) 
      stop("Host D matrix has negative eigenvalues. Rerun with correction=\"lingoes\" or correction=\"cailliez\"")
    sum.host.values.sq <- sum(host.pc$values[, 1]^2)
    host.vectors <- host.pc$vectors
  }else {
    sum.host.values.sq <- sum(host.pc$values[, 2]^2)
    host.vectors <- host.pc$vectors.cor
  }
  if (!nullmodel %in% c(1,2)) {
    stop("Pick a null model among \"1\" or \"2\".")
  }
  n.host <- nrow(host.D)
  para.D <- as.matrix(para.D)
  para.pc <- pcoa_test(para.D, correction = correction)
  if (para.pc$correction[2] == 1) {
    if (min(para.pc$values[, 2]) < -epsilon) 
      stop("Parasite D matrix has negative eigenvalues. Rerun with correction=\"lingoes\" or correction=\"cailliez\"")
    sum.para.values.sq <- sum(para.pc$values[, 1]^2)
    para.vectors <- para.pc$vectors
  }else {
    sum.para.values.sq <- sum(para.pc$values[, 2]^2)
    para.vectors <- para.pc$vectors.cor
  }
  n.para <- nrow(para.D)
  if (!silent) 
    cat("n.hosts =", n.host, ", n.parasites =", n.para, "\n")
  a <- system.time({
    tracemax <- max(sum.host.values.sq, sum.para.values.sq)
    if (n.host == n.para) {
      if (!silent) 
        cat("The function cannot check if matrix HP has been entered in the right way.", 
            "\n")
      if (!silent) 
        cat("It will assume that the rows of HP are the hosts.", 
            "\n")
    }else {
      temp <- dim(HP)
      if (temp[1] == n.host) {
        if (temp[2] != n.para) 
          stop("Matrices host.D, para.D and HP not comformable")
      }else if (temp[2] == n.host) {
        if (temp[1] != n.para) 
          stop("Matrices host.D, para.D and HP not comformable")
        HP <- t(HP)
        if (!silent) 
          cat("Matrix HP has been transposed for comformity with host.D and para.D.", 
              "\n")
      }else {
        stop("Matrices host.D, para.D and HP not comformable")
      }
    }
    p.per.h <- apply(HP, 1, sum)
    h.per.p <- apply(HP, 2, sum)
    mat.4 <- t(host.vectors) %*% HP %*% para.vectors
    global <- sum(mat.4^2)
    if (nperm > 0) {
      set.seed(seed)
      nGT <- 1
      global.perm <- NA
      for (i in 1:nperm) {
        if (nullmodel==1) {HP.perm <- apply(HP, 2, sample)}
        if (nullmodel==2) {
          HP.perm <- HP[sample(1:nrow(HP)),]
        }
        mat.4.perm <- t(host.vectors) %*% HP.perm %*% para.vectors
        global.perm <- c(global.perm, sum(mat.4.perm^2))
        if (global.perm[i + 1] >= global) 
          nGT <- nGT + 1
      }
      global.perm <- global.perm[-1]
      p.global <- nGT/(nperm + 1)
    }else {
      p.global <- NA
    }
    if (test.links) {
      list.hp <- which(t(cbind(HP, rep(0, n.host))) > 0)
      HP.list <- cbind((list.hp%/%(n.para + 1)) + 1, list.hp%%(n.para + 1))
      colnames(HP.list) <- c("Host", "Parasite")
      n.links <- length(list.hp)
      stat1 <- NA
      stat2 <- NA
      p.stat1 <- NA
      p.stat2 <- NA
      for (k in 1:n.links) {
        HP.k <- HP
        HP.k[HP.list[k, 1], HP.list[k, 2]] <- 0
        mat.4.k <- t(host.vectors) %*% HP.k %*% para.vectors
        trace.k <- sum(mat.4.k^2)
        stat1 <- c(stat1, (global - trace.k))
        den <- tracemax - global
        if (den > epsilon) {
          stat2 <- c(stat2, stat1[k + 1]/den)
        }else {
          stat2 <- c(stat2, NA)
        }
        if (nperm > 0) {
          set.seed(seed)
          nGT1 <- 1
          nGT2 <- 1
          nperm2 <- nperm
          for (i in 1:nperm) {
            HP.k.perm <- apply(HP.k, 2, sample)
            mat.4.k.perm <- t(host.vectors) %*% HP.k.perm %*% 
              para.vectors
            trace.k.perm <- sum(mat.4.k.perm^2)
            stat1.perm <- global.perm[i] - trace.k.perm
            if (stat1.perm >= stat1[k + 1]) 
              nGT1 <- nGT1 + 1
            if (!is.na(stat2[k + 1])) {
              den <- tracemax - global.perm[i]
              if (den > epsilon) {
                stat2.perm <- stat1.perm/den
                if (stat2.perm >= stat2[k + 1]) 
                  nGT2 <- nGT2 + 1
              }else {
                nperm2 <- nperm2 - 1
              }
            }
          }
          p.stat1 <- c(p.stat1, nGT1/(nperm + 1))
          if (!is.na(stat2[k + 1])) {
            p.stat2 <- c(p.stat2, nGT2/(nperm2 + 1))
          }else {
            p.stat2 <- c(p.stat2, NA)
          }
        }else {
          p.stat1 <- c(p.stat1, NA)
          p.stat2 <- c(p.stat2, NA)
        }
      }
      link.table <- cbind(HP.list, stat1[-1], p.stat1[-1], 
                          stat2[-1], p.stat2[-1])
      colnames(link.table) = c("Host", "Parasite", "F1.stat", 
                               "p.F1", "F2.stat", "p.F2")
      out <- list(ParaFitGlobal = global, p.global = p.global, 
                  link.table = link.table, para.per.host = p.per.h, 
                  host.per.para = h.per.p, nperm = nperm)
    }else {
      if (!silent) 
        cat("Rerun the program with option 'test.links=TRUE' to test the individual H-P links", 
            "\n")
      out <- list(ParaFitGlobal = global, p.global = p.global, 
                  para.per.host = p.per.h, host.per.para = h.per.p, 
                  nperm = nperm)
    }
  })
  a[3] <- sprintf("%2f", a[3])
  if (!silent) 
    cat("Computation time =", a[3], " sec", "\n")
  class(out) <- "parafit"
  out
}


pcoa_test <-  function (D, correction = "none", rn = NULL) {
  centre <- function(D, n) {
    One <- matrix(1, n, n)
    mat <- diag(n) - One/n
    mat.cen <- mat %*% D %*% mat
  }
  bstick.def <- function(n, tot.var = 1, ...) {
    res <- rev(cumsum(tot.var/n:1)/n)
    names(res) <- paste("Stick", seq(len = n), sep = "")
    return(res)
  }
  D <- as.matrix(D)
  n <- nrow(D)
  epsilon <- sqrt(.Machine$double.eps)
  if (length(rn) != 0) {
    names <- rn
  } else {
    names <- rownames(D)
  }
  CORRECTIONS <- c("none", "lingoes", "cailliez")
  correct <- pmatch(correction, CORRECTIONS)
  if (is.na(correct)) 
    stop("Invalid correction method")
  delta1 <- centre((-0.5 * D^2), n)
  trace <- sum(diag(delta1))
  D.eig <- eigen(delta1)
  min.eig <- min(D.eig$values)
  zero.eig <- which(abs(D.eig$values) < epsilon)
  D.eig$values[zero.eig] <- 0
  
  if (min.eig > -epsilon) {
    correct <- 1
    eig <- D.eig$values
    k <- length(which(eig > epsilon))
    rel.eig <- eig[1:k]/trace
    cum.eig <- cumsum(rel.eig)
    vectors <- sweep(D.eig$vectors[, 1:k, drop=F], 2, sqrt(eig[1:k]), 
                     FUN = "*")
    bs <- bstick.def(k)
    cum.bs <- cumsum(bs)
    res <- data.frame(eig[1:k], rel.eig, bs, cum.eig, cum.bs)
    colnames(res) <- c("Eigenvalues", "Relative_eig", "Broken_stick", 
                       "Cumul_eig", "Cumul_br_stick")
    rownames(res) <- 1:nrow(res)
    rownames(vectors) <- names
    colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                                  prefix = "Axis.")
    note <- paste("There were no negative eigenvalues. No correction was applied")
    out <- (list(correction = c(correction, correct), note = note, 
                 values = res, vectors = vectors, trace = trace))
  }else {
    k <- n
    eig <- D.eig$values
    rel.eig <- eig/trace
    rel.eig.cor <- (eig - min.eig)/(trace - (n - 1) * min.eig)
    if (length(zero.eig)) 
      rel.eig.cor <- c(rel.eig.cor[-zero.eig[1]], 0)
    cum.eig.cor <- cumsum(rel.eig.cor)
    k2 <- length(which(eig > epsilon))
    k3 <- length(which(rel.eig.cor > epsilon))
    vectors <- sweep(D.eig$vectors[, 1:k2, drop=F], 2, sqrt(eig[1:k2]), 
                     FUN = "*")
    if ((correct == 2) | (correct == 3)) {
      if (correct == 2) {
        c1 <- -min.eig
        note <- paste("Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 -", 
                      c1, ", except diagonal elements")
        D <- -0.5 * (D^2 + 2 * c1)
      }else if (correct == 3) {
        delta2 <- centre((-0.5 * D), n)
        upper <- cbind(matrix(0, n, n), 2 * delta1)
        lower <- cbind(-diag(n), -4 * delta2)
        sp.matrix <- rbind(upper, lower)
        c2 <- max(Re(eigen(sp.matrix, symmetric = FALSE, 
                           only.values = TRUE)$values))
        note <- paste("Cailliez correction applied to negative eigenvalues: D' = -0.5*(D +", 
                      c2, ")^2, except diagonal elements")
        D <- -0.5 * (D + c2)^2
      }
      diag(D) <- 0
      mat.cor <- centre(D, n)
      toto.cor <- eigen(mat.cor)
      trace.cor <- sum(diag(mat.cor))
      min.eig.cor <- min(Re(toto.cor$values)) # min(toto.cor$values) # problem with imaginary number close to 0
      zero.eig.cor <- which((Re(toto.cor$values) < epsilon) & 
                              (Re(toto.cor$values) > -epsilon)) # idem
      toto.cor$values[zero.eig.cor] <- 0
      toto.cor$values <- Re(toto.cor$values)
      toto.cor$vectors <- Re(toto.cor$vectors) # idem
      
      if (min.eig.cor > -epsilon) {
        eig.cor <- toto.cor$values
        rel.eig.cor <- eig.cor[1:k]/trace.cor
        cum.eig.cor <- cumsum(rel.eig.cor)
        k2 <- length(which(eig.cor > epsilon))
        vectors.cor <- sweep(toto.cor$vectors[, 1:k2,drop=F], 
                             2, sqrt(eig.cor[1:k2]), FUN = "*")
        rownames(vectors.cor) <- names
        colnames(vectors.cor) <- colnames(vectors.cor, 
                                          do.NULL = FALSE, prefix = "Axis.")
        bs <- bstick.def(k2)
        bs <- c(bs, rep(0, (k - k2)))
        cum.bs <- cumsum(bs)
      }else {
        if (correct == 2) 
          cat("Problem! Negative eigenvalues are still present after Lingoes", 
              "\n")
        if (correct == 3) 
          cat("Problem! Negative eigenvalues are still present after Cailliez", 
              "\n")
        rel.eig.cor <- cum.eig.cor <- bs <- cum.bs <- rep(NA, 
                                                          n)
        vectors.cor <- matrix(NA, n, 2)
        rownames(vectors.cor) <- names
        colnames(vectors.cor) <- colnames(vectors.cor, 
                                          do.NULL = FALSE, prefix = "Axis.")
      }
      res <- data.frame(eig[1:k], eig.cor[1:k], rel.eig.cor, 
                        bs, cum.eig.cor, cum.bs)
      colnames(res) <- c("Eigenvalues", "Corr_eig", "Rel_corr_eig", 
                         "Broken_stick", "Cum_corr_eig", "Cum_br_stick")
      rownames(res) <- 1:nrow(res)
      rownames(vectors) <- names
      colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                                    prefix = "Axis.")
      out <- (list(correction = c(correction, correct), 
                   note = note, values = res, vectors = vectors, 
                   trace = trace, vectors.cor = vectors.cor, trace.cor = trace.cor))
    }else {
      note <- "No correction was applied to the negative eigenvalues"
      bs <- bstick.def(k3)
      bs <- c(bs, rep(0, (k - k3)))
      cum.bs <- cumsum(bs)
      res <- data.frame(eig[1:k], rel.eig, rel.eig.cor, 
                        bs, cum.eig.cor, cum.bs)
      colnames(res) <- c("Eigenvalues", "Relative_eig", 
                         "Rel_corr_eig", "Broken_stick", "Cum_corr_eig", 
                         "Cumul_br_stick")
      rownames(res) <- 1:nrow(res)
      rownames(vectors) <- names
      colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                                    prefix = "Axis.")
      out <- (list(correction = c(correction, correct), 
                   note = note, values = res, vectors = vectors, 
                   trace = trace))
    }
  }
  class(out) <- "pcoa"
  out
}


########  functions from the packages PACo ####


add_pcoord_test <- function (D) {
  HP_bin <- which(D$HP > 0, arr.ind = TRUE)
  H_PCo <- coordpcoa_test(D$H, correction = "cailliez")$vectors
  P_PCo <- coordpcoa_test(D$P, correction = "cailliez")$vectors
  D$H_PCo <- H_PCo[HP_bin[, 1,drop=F], ,drop=F]
  D$P_PCo <- P_PCo[HP_bin[, 2,drop=F], ,drop=F]
  return(D)
}

coordpcoa_test <- function (D, correction = "none", rn = NULL) {
  centre <- function(D, n) {
    One <- matrix(1, n, n)
    mat <- diag(n) - One/n
    mat.cen <- mat %*% D %*% mat
  }
  bstick.def <- function(n, tot.var = 1, ...) {
    res <- rev(cumsum(tot.var/n:1)/n)
    names(res) <- paste("Stick", seq(len = n), sep = "")
    return(res)
  }
  D <- as.matrix(D)
  n <- nrow(D)
  epsilon <- sqrt(.Machine$double.eps)
  if (length(rn) != 0) {
    names <- rn
  }else {
    names <- rownames(D)
  }
  CORRECTIONS <- c("none", "lingoes", "cailliez")
  correct <- pmatch(correction, CORRECTIONS)
  if (is.na(correct)) 
    stop("Invalid correction method")
  delta1 <- centre((-0.5 * D^2), n)
  trace <- sum(diag(delta1))
  D.eig <- eigen(delta1)
  D.eig$values <- as.numeric(zapsmall(vegan::eigenvals(D.eig)))
  min.eig <- min(D.eig$values)
  zero.eig <- which(abs(D.eig$values) < epsilon)
  D.eig$values[zero.eig] <- 0
  if (min.eig > -epsilon) {
    correct <- 1
    eig <- D.eig$values
    k <- length(which(eig > epsilon))
    rel.eig <- eig[1:k]/trace
    cum.eig <- cumsum(rel.eig)
    vectors <- sweep(D.eig$vectors[, 1:k,drop=F], 2, sqrt(eig[1:k,drop=F]), 
                     FUN = "*")
    bs <- bstick.def(k)
    cum.bs <- cumsum(bs)
    res <- data.frame(eig[1:k], rel.eig, bs, cum.eig, cum.bs)
    colnames(res) <- c("Eigenvalues", "Relative_eig", "Broken_stick", 
                       "Cumul_eig", "Cumul_br_stick")
    rownames(res) <- 1:nrow(res)
    rownames(vectors) <- names
    colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                                  prefix = "Axis.")
    note <- paste("There were no negative eigenvalues. No correction was applied")
    out <- (list(correction = c(correction, correct), note = note, 
                 values = res, vectors = vectors, trace = trace))
  }else {
    k <- n
    eig <- D.eig$values
    rel.eig <- eig/trace
    rel.eig.cor <- (eig - min.eig)/(trace - (n - 1) * min.eig)
    rel.eig.cor = c(rel.eig.cor[1:(zero.eig[1] - 1)], rel.eig.cor[(zero.eig[1] + 
                                                                     1):n], 0)
    cum.eig.cor <- cumsum(rel.eig.cor)
    k2 <- length(which(eig > epsilon))
    k3 <- length(which(rel.eig.cor > epsilon))
    vectors <- sweep(D.eig$vectors[, 1:k2,drop=F], 2, sqrt(eig[1:k2,drop=F]), 
                     FUN = "*")
    if ((correct == 2) | (correct == 3)) {
      if (correct == 2) {
        c1 <- -min.eig
        note <- paste("Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 -", 
                      c1, ", except diagonal elements")
        D <- -0.5 * (D^2 + 2 * c1)
      }else if (correct == 3) {
        delta2 <- centre((-0.5 * D), n)
        upper <- cbind(matrix(0, n, n), 2 * delta1)
        lower <- cbind(-diag(n), -4 * delta2)
        sp.matrix <- rbind(upper, lower)
        c2 <- max(Re(eigen(sp.matrix, symmetric = FALSE, 
                           only.values = TRUE)$values))
        note <- paste("Cailliez correction applied to negative eigenvalues: D' = -0.5*(D +", 
                      c2, ")^2, except diagonal elements")
        D <- -0.5 * (D + c2)^2
      }
      diag(D) <- 0
      mat.cor <- centre(D, n)
      toto.cor <- eigen(mat.cor)
      toto.cor$values <- as.numeric(zapsmall(vegan::eigenvals(toto.cor)))
      trace.cor <- sum(diag(mat.cor))
      min.eig.cor <- min(toto.cor$values)
      zero.eig.cor <- which((toto.cor$values < epsilon) & 
                              (toto.cor$values > -epsilon))
      toto.cor$values[zero.eig.cor] <- 0
      if (min.eig.cor > -epsilon) {
        eig.cor <- toto.cor$values
        rel.eig.cor <- eig.cor[1:k]/trace.cor
        cum.eig.cor <- cumsum(rel.eig.cor)
        k2 <- length(which(eig.cor > epsilon))
        vectors.cor <- sweep(toto.cor$vectors[, 1:k2,drop=F], 
                             2, sqrt(eig.cor[1:k2,drop=F]), FUN = "*")
        rownames(vectors.cor) <- names
        colnames(vectors.cor) <- colnames(vectors.cor, 
                                          do.NULL = FALSE, prefix = "Axis.")
        bs <- bstick.def(k2)
        bs <- c(bs, rep(0, (k - k2)))
        cum.bs <- cumsum(bs)
      }else {
        if (correct == 2) 
          cat("Problem! Negative eigenvalues are still present after Lingoes", 
              "\n")
        if (correct == 3) 
          cat("Problem! Negative eigenvalues are still present after Cailliez", 
              "\n")
        rel.eig.cor <- cum.eig.cor <- bs <- cum.bs <- rep(NA, 
                                                          n)
        vectors.cor <- matrix(NA, n, 2)
        rownames(vectors.cor) <- names
        colnames(vectors.cor) <- colnames(vectors.cor, 
                                          do.NULL = FALSE, prefix = "Axis.")
      }
      res <- data.frame(eig[1:k], eig.cor[1:k], rel.eig.cor, 
                        bs, cum.eig.cor, cum.bs)
      colnames(res) <- c("Eigenvalues", "Corr_eig", "Rel_corr_eig", 
                         "Broken_stick", "Cum_corr_eig", "Cum_br_stick")
      rownames(res) <- 1:nrow(res)
      rownames(vectors) <- names
      colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                                    prefix = "Axis.")
      out <- (list(correction = c(correction, correct), 
                   note = note, values = res, vectors = vectors, 
                   trace = trace, vectors.cor = vectors.cor, trace.cor = trace.cor))
    }else {
      note <- "No correction was applied to the negative eigenvalues"
      bs <- bstick.def(k3)
      bs <- c(bs, rep(0, (k - k3)))
      cum.bs <- cumsum(bs)
      res <- data.frame(eig[1:k], rel.eig, rel.eig.cor, 
                        bs, cum.eig.cor, cum.bs)
      colnames(res) <- c("Eigenvalues", "Relative_eig", 
                         "Rel_corr_eig", "Broken_stick", "Cum_corr_eig", 
                         "Cumul_br_stick")
      rownames(res) <- 1:nrow(res)
      rownames(vectors) <- names
      colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                                    prefix = "Axis.")
      out <- (list(correction = c(correction, correct), 
                   note = note, values = res, vectors = vectors, 
                   trace = trace))
    }
  }
  class(out) <- "pcoa"
  out
}


PACo_test <- function (D, nperm = 1000, seed = NA, nullmodel = 1) {
  if (!nullmodel %in% c(1,2)) {
    stop("Pick a null model among \"1\" or \"2\".")
  }
  if (!("H_PCo" %in% names(D))) 
    D <- add_pcoord_test(D)
  proc <- vegan::procrustes(X = D$H_PCo, Y = D$P_PCo)
  Nlinks <- sum(D$HP)
  m2ss <- proc$ss
  pvalue <- 0
  if (!is.na(seed)) 
    set.seed(seed)
  
  if (nullmodel==1){
    null_model <- vegan::nullmodel(D$HP, "r0")
    randomised_matrices <- stats::simulate(null_model, nsim = nperm)
  }
  
  for (n in c(1:nperm)) {
    if (nullmodel==1) {
      permuted_HP <- randomised_matrices[, , n]
      permuted_HP <- permuted_HP[rownames(D$HP), colnames(D$HP)]
    } else { # null model 2
      permuted_HP <- D$HP[sample(1:nrow(D$HP)),]
      rownames(permuted_HP) <- rownames(D$HP)
    }
    perm_D <- list(H = D$H, P = D$P, HP = permuted_HP)
    perm_paco <- add_pcoord_test(perm_D)
    perm_proc_ss <- vegan::procrustes(X = perm_paco$H_PCo, 
                                      Y = perm_paco$P_PCo)$ss
    
    if (is.na(perm_proc_ss)){
      nperm <- nperm - 1
    } else if (perm_proc_ss <= m2ss) {
      pvalue <- pvalue + 1
    }
  }
  pvalue <- pvalue/nperm
  D$proc <- proc
  D$gof <- list(p = pvalue, ss = m2ss, n = nperm)
  D$nullmodel <- nullmodel
  return(D)
}




