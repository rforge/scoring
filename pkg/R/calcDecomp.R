calcDecomp <- function(f, d, bin, wt, nbin = max(bin), scale=FALSE, ...){
  ## compute weighted brier decomposition
  ## nbin is total number of bins across all systems
  dotdotdot <- list(...)
  useifps <- "ifpid" %in% names(dotdotdot)
  if(useifps){
    ifpid <- dotdotdot$ifpid
  }

  ## what bins are missing for this system?
  ## add phantom observations for these bins with weights of 0
  misbins <- which(!(1:nbin %in% bin))
  if(length(misbins) > 0){
    bin <- c(bin, misbins)
    tmpf <- matrix(1/ncol(f), length(misbins), ncol(f))
    colnames(tmpf) <- colnames(f)
    f <- rbind(f, tmpf)
    tmpd <- matrix(0, length(misbins), ncol(f))
    colnames(tmpd) <- colnames(d)
    d <- rbind(d, tmpd)
    wt <- c(wt, rep(0, length(misbins)))
    if(useifps) ifpid <- c(ifpid, rep(ifpid[1], length(misbins)))
  }

  binw <- tapply(wt, bin, sum) # w* from paper
  dw <- apply(d, 2, function(x) x*wt)
  dbar <- apply(dw, 2, sum)
  bindw <- apply(dw, 2, function(x) tapply(x, bin, sum)/binw) # dbar_km from paper
  binf <- apply(f, 2, function(x) tapply(x, bin, head, 1))
  
  unc <- sum(dbar * (1 - dbar))

  if(useifps){
    ## assume weights sum to 1, each ifpid weighted equally
    bindbar <- bindw[bin,]
    qdiscrim <- apply(bindbar, 1, function(x) sum((x - dbar)^2, na.rm = TRUE))

    qwt <- tapply(wt, ifpid, sum)
    qwt <- qwt[match(ifpid, as.numeric(names(qwt)))]
    qdiscrim <- tapply(wt * qdiscrim / qwt, ifpid, sum)
    discrim <- mean(qdiscrim)
  } else {
    ## weights don't necessarily sum to 1
    qdiscrim <- NULL
    discrim <- sum(binw * apply(bindw, 1, function(x) sum((x - dbar)^2, na.rm = TRUE)))
  }

  if(useifps){
    ## calibration
    qcalib <- apply((f - bindbar)^2, 1, sum, na.rm = TRUE)
    qcalib <- tapply(wt * qcalib / qwt, ifpid, sum)
    calib <- mean(qcalib)
  } else {
    ## weights don't necessarily sum to 1
    qcalib <- NULL
    calib <- sum(binw * apply((binf - bindw)^2, 1, sum, na.rm = TRUE))
  }

  ## usual mmde calculation:
  mmde <- sum(wt * apply((f - d)^2, 1, sum))
  mmdecomp <- unc - discrim + calib

  if(abs(mmde - mmdecomp) > 1e-5){
    print(c(mmde, mmdecomp))
    stop("Problem with decomposition (does not sum to MMDE).")
  }

  ## yates decomposition
  fbar <- apply(f, 2, function(x) sum(wt * x))
  fcent <- t(apply(f, 1, function(x) (x - fbar)))
  dcent <- t(apply(d, 1, function(x) (x - dbar)))
  fvar <- apply(fcent, 2, function(x) sum(wt * x^2))
  covfd <- fcent * dcent
  covfd <- apply(covfd, 2, function(x) sum(wt * x))
  ##covfd2 <- rep(NA, ncol(f))
  ##for(j in 1:ncol(f)){
  ##  covfd2[j] <- sum(wt * (f[,j] - fbar[j]) * (d[,j] - dbar[j]))
  ##}
  disc2 <- t(apply(bindw, 1, function(x) (x - dbar)^2))
  disc2 <- apply(disc2, 2, function(x) sum(binw * x, na.rm = TRUE))
  calib2 <- sum((fbar - dbar)^2)

  dna <- d; dna[d == 0] <- NA
  f1 <- f * dna
  dna <- d; dna[d == 1] <- NA
  f0 <- f * (1 - dna)
  fbar1 <- apply(f1, 2, function(x) sum(wt * x, na.rm = TRUE) /
                                    sum(wt * !is.na(x)))
  ## in case no d=1 for an alternative
  fbar1[is.na(fbar1)] <- 0
  fbar0 <- apply(f0, 2, function(x) sum(wt * x, na.rm = TRUE) /
                                    sum(wt * !is.na(x)))
  ## in case no d=0 for an alternative
  fbar0[is.na(fbar0)] <- 0
  deltaf <- sum(fvar - (fbar1 - fbar0)^2 * dbar * (1 - dbar))

  ## check per yates 1982, eq (10)
  calibalt <- sum(fvar + (fbar - dbar)^2 - 2*covfd + disc2)
  if(abs(calib - calibalt) > 1e-5){
    stop("Problem with Yates decomposition.")
  }
  
  res <- c(discrim, calib, deltaf, calib2, sum(covfd), unc)

  if(scale){
    ubnd <- (ncol(f1) - 1)/ncol(f1)
    compmin <- rep(0,6)
    compmax <- c(ubnd, 2, sum(fvar), 1, ubnd, ubnd)

    res <- (res - compmin)/compmax
    ## reverse cal/var metrics so bigger is better
    res[2:4] <- 1 - res[2:4]
  }
  if(length(res) > 6) stop("Problem with decomposition summary.")

  list(res = res, qcalib = qcalib, qdiscrim = qdiscrim)
}

