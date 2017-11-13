## scoring helper functions

## add phantom alternatives, change ordinal questions
## into binary, cumulative probs
phOifp <- function(dat, ngrp, fcols, qtype){
  fdat <- vector("list", ngrp)
  ddat <- vector("list", ngrp)
  wtdat <- vector("list", ngrp)

  qinfo <- !is.null(qtype)

  if(qinfo) ifpdat <- vector("list", ngrp)
  for(i in 1:ngrp){
    fdat[[i]] <- dat[[i]][,fcols]

    if(qinfo){
      ifpid <- dat[[i]]$qid
      ifpmatch <- match(ifpid, qtype$qid)
      ifpdat[[i]] <- cbind.data.frame(ifpid,
                                      squo = qtype$squo[ifpmatch],
                                      ord = qtype$ord[ifpmatch])
    }

    ddat[[i]] <- matrix(0, nrow(fdat[[i]]), ncol(fdat[[i]]))
    ddat[[i]][cbind(1:nrow(fdat[[i]]),
                    as.numeric(dat[[i]][["outcome"]]))] <- 1

    wtdat[[i]] <- dat[[i]][["wt"]]
    
    if(any(qtype$ord == 1)){
      orddat <- expandOrd(fdat[[i]], ddat[[i]], ifpdat[[i]],
                          wtdat[[i]])

      fdat[[i]] <- orddat$fdat
      ddat[[i]] <- orddat$ddat
      ifpdat[[i]] <- orddat$ifpdat
      wtdat[[i]] <- orddat$wtdat
    }
    
    fdat[[i]][is.na(fdat[[i]])] <- 0
    ddat[[i]][is.na(ddat[[i]])] <- 0
  }

  list(fdat = fdat, ddat = ddat, ifpdat = ifpdat, wtdat = wtdat)
}
  

## expand out ordinal ifps into cumulatives
expandOrd <- function(fdat, ddat, ifpdat, wtdat){
  oifps <- which(ifpdat$ord == 1)

  cumfs <- t(apply(fdat[oifps,], 1, cumsum))
  cumds <- t(apply(ddat[oifps,], 1, cumsum))
  nalts <- apply(!is.na(fdat[oifps,]), 1, sum)

  M <- ncol(fdat)
  newifpid <- rep(ifpdat$ifpid[oifps], (nalts-1))
  neword <- rep(ifpdat$ord[oifps], (nalts-1))
  newsquo <- rep(ifpdat$squo[oifps], (nalts-1))
  newwts <- rep(wtdat[oifps] / (nalts-1), (nalts-1))
  
  for(i in 1:length(oifps)){
    cumfs[i, nalts[i]:M] <- NA
    cumds[i, nalts[i]:M] <- NA
  }

  cumfvec <- as.numeric(t(cumfs))
  cumfvec <- as.numeric(na.omit(cumfvec))
  cumdvec <- as.numeric(t(cumds))
  cumdvec <- as.numeric(na.omit(cumdvec))

  ## first two cols are cumulative probs, the remaining are phantoms
  newf <- cbind(cumfvec, 1-cumfvec, matrix(NA, length(cumfvec), M-2))
  colnames(newf) <- names(fdat)
  newd <- cbind(cumdvec, 1-cumdvec, matrix(NA, length(cumdvec), M-2))
  colnames(newd) <- names(ddat)
  
  fdat <- rbind(fdat[-oifps,], newf)
  ddat <- rbind(ddat[-oifps,], newd)

  ifpdat <- rbind(ifpdat[-oifps,], cbind.data.frame(ifpid=newifpid,
                                                    ord=neword,
                                                    squo=newsquo))

  wtdat <- c(wtdat[-oifps], newwts)
  
  list(fdat = fdat, ddat = ddat, ifpdat = ifpdat, wtdat = wtdat)
}


setBins <- function(fdat, bin){
  grp <- rep(1:length(fdat), sapply(fdat, nrow))
  fullf <- do.call("rbind", fdat)
  if(bin){
    fullf <- round(fullf, 1)
    binfs <- t(apply(fullf, 1, function(x){
      if(sum(x, na.rm=TRUE) != 1){
        mv <- which(x == min(x, na.rm = TRUE))
        if(length(mv) > 1) mv <- sample(mv, 1)
        x[mv] <- 1 - sum(x[-mv], na.rm = TRUE)
      }
      x
    }))
  } else {
    binfs <- fullf
  }
  
  ## separate out fs by system
  binfs <- split.data.frame(binfs, grp)

  binfs
}


## shuffle forecast alternatives separately for
##  status quo ifps (all shuffled identically)
##  ordinal ifps (all shuffled identically)
##  unordered ifps (each shuffled differently)
## Note if !shuff, then we just organize data for
## calcDecomp()
shuffle <- function(fdat, ddat, ifpdat, shuff = TRUE){
  nalt <- ncol(fdat[[1]])
  nsystem <- length(fdat)

  if(shuff){
    sqshuf <- sample(1:nalt, nalt)
    oshuf <- sample(1:nalt, nalt)
  } else {
    sqshuf <- 1:nalt
    oshuf <- 1:nalt
  }

  for(i in 1:nsystem){
    if(any(ifpdat[[i]][['squo']] == 1 & ifpdat[[i]][['ord']] ==1)) stop("IFPs cannot be coded as both 'status quo' and 'ordinal'.")

    unifp <- unique(ifpdat[[i]][['ifpid']][ifpdat[[i]][['squo']] != 1 & ifpdat[[i]][['ord']] != 1])
    
    ## a single shuffle for squo questions, and for ordinal questions
    fdat[[i]][ifpdat$squo == 1,] <- fdat[[i]][ifpdat$squo == 1, sqshuf]
    fdat[[i]][ifpdat$ord == 1,] <- fdat[[i]][ifpdat$ord == 1, oshuf]
    ddat[[i]][ifpdat$squo == 1,] <- ddat[[i]][ifpdat$squo == 1, sqshuf]
    ddat[[i]][ifpdat$ord == 1,] <- ddat[[i]][ifpdat$ord == 1, oshuf]
  }

  tmpshuf <- 1:nalt
  for(j in 1:length(unifp)){
    if(shuff) tmpshuf <- sample(1:nalt, nalt)
    for(i in 1:nsystem){
      tmprows <- ifpdat[[i]][['ifpid']] == unifp[j]
      fdat[[i]][tmprows,1:nalt] <- fdat[[i]][tmprows,tmpshuf]
      ddat[[i]][tmprows,1:nalt] <- ddat[[i]][tmprows,tmpshuf]
    }
  }

  ## define new bins based on shuffled rows (without
  ## rounding/transforming any forecasts)
  fullf <- do.call("rbind", fdat)
  newbins <- paste0(fullf[,1], fullf[,2])
  if(ncol(fullf) > 2){
    for(i in 3:ncol(fullf)){
      newbins <- paste0(newbins, fullf[,i])
    }
  }
  newbins <- as.numeric(factor(newbins))
  nbin <- max(newbins)
  newbins <- split(newbins, rep(1:nsystem, each=nrow(fullf)/nsystem))
  
  list(fdat=fdat, ddat=ddat, newbins=newbins, nbin=nbin)
}
