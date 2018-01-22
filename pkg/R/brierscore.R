brierscore <- function(object, data, group=NULL, decomp=FALSE, bounds=NULL, reverse=FALSE, wt=NULL, decompControl=list()){

    ## Machinery to get group variable, if specified
    ## FIXME Use | within a formula for groups?
    if(missing(data)) data <- environment(object)
    ## So bounds of $rawscores match the decomposition:
    ##if(decomp & is.null(bounds)) bounds <- c(0,2)

    ## ensure wt + decompControl$wt are the same
    ## (wt arg only exists in case user wants a weighted
    ##  brier score without doing the decomposition)
    if(!is.null(wt) & ("wt" %in% names(decompControl))){
        if(!(all.equal(wt, decompControl$wt))){
            stop("Different weights supplied to wt and decompControl arguments.")
        }
    }
    if(!is.null(wt) & !("wt" %in% names(decompControl))){
        decompControl <- c(decompControl, list(wt=wt))
    }
    if(is.null(wt) & ("wt" %in% names(decompControl))){
        wt <- decompControl$wt
    }

    if("qtype" %in% names(decompControl)){
        qtype <- decompControl$qtype
        ordrows <- as.logical(qtype$ord[match(data$question.id,
                                              qtype$qid)])
    } else {
        ordrows <- FALSE
    }
        
    args <- list(object=object, param=2, data=data, bounds=bounds,
                 reverse=reverse, ordered = ordrows, decomp=decomp,
                 decompControl=decompControl)

    if(!is.null(group)){
        mf <- match.call()
        m <- as.character(mf[match("group", names(mf), 0L)])
        group <- data[[m]]
        args <- c(args, list(group=group))
    }

    ## Get Brier scores + decomp using calcscore
    scores <- do.call("calcscore", args)

    ## If they supply group but not decomp, at least give them
    ## avg Briers
    if(!decomp & !is.null(group)){
        if(is.null(wt)){
            mnbrier <- tapply(scores, group, mean)
        } else {
            mnbrier <- tapply(scores*wt, group, sum)
        }
        scores <- list(rawscores=scores, brieravg=mnbrier)
    }

    scores
}

bdecomp <- function(forecast, outcome,
                    group = rep(1, length(outcome)),
                    wt = NULL, qid = 1:nrow(forecast), bin = TRUE,
                    qtype = NULL, scale = FALSE,
                    roundto = .1, binstyle = 1,
                    resamples = 0){

    origqtype <- qtype
    if(is.null(wt)){
        wt <- 1/table(group)
        wttrans <- match(as.character(group), names(wt))
        wt <- as.numeric(wt[wttrans])
    }

    ## ensure forecasts sum to 1
    fsums <- data.frame(t(apply(forecast, 1, function(x) x/sum(x, na.rm = TRUE))))
    names(fsums) <- names(forecast)
    forecast <- fsums

    dat <- cbind.data.frame(forecast, outcome = outcome,
                            group = group, wt = wt)
    if(is.null(qtype)){
        uqs <- unique(qid)
        qtype <- cbind.data.frame(qid = uqs,
                                  ord = rep(0, length(uqs)),
                                  squo = rep(0, length(uqs)))

    }

    dat <- cbind.data.frame(dat, qid = qid)
  
    dat <- split.data.frame(dat, group)
  
    ngrp <- length(dat)
    grpnames <- names(dat)
  
    ## compute original (unrounded) mmde
    if(bin){
        origmmde <- bdecomp(forecast, outcome, group, wt, qid,
                            bin=FALSE, qtype=origqtype, scale=scale)
    }

    ## add phantoms, recode oIFPs
    newdat <- phOifp(dat, ngrp, 1:ncol(forecast), qtype)

    n <- sapply(newdat$fdat, nrow)
    J <- sapply(newdat$ifpdat, function(x) length(unique(x$ifpid)))
    if(length(unique(n)) > 1) stop("Different groups have differing numbers of forecasts. Group comparisons will be inaccurate.")  

    binfs <- setBins(newdat$fdat, bin, roundto, binstyle)

    reslist <- qcalib <- qdiscrim <- vector("list", ngrp)
    niter <- ifelse(resamples < 1, 1, resamples)
    for(i in 1:ngrp){
      reslist[[i]] <- matrix(NA, niter, (6 + !bin))
      qcalib[[i]] <- matrix(NA, niter, J)
      qdiscrim[[i]] <- matrix(NA, niter, J)
    }


    for(i in 1:niter){
        shuffdat <- shuffle(binfs, newdat$ddat, newdat$ifpdat,
                            shuff = resamples > 1)
            
        for(j in 1:ngrp){
            tmpf <- shuffdat$fdat[[j]]
            tmpd <- shuffdat$ddat[[j]]
            tmpbin <- shuffdat$newbins[[j]]
            nbin <- shuffdat$nbin
            tmpwt <- newdat$wtdat[[j]]

            if(bin){
                decomp <- calcDecomp(tmpf, tmpd, tmpbin, tmpwt, nbin, scale = scale, ifpid = newdat$ifpdat[[j]]$ifpid)

                reslist[[j]][i,] <- decomp$res
                qcalib[[j]][i,] <- decomp$qcalib
                qdiscrim[[j]][i,] <- decomp$qdiscrim
            } else {
                rowbriers <- apply((tmpf - tmpd)^2, 1, sum)
                mde <- tapply(tmpwt * rowbriers, newdat$ifpdat[[j]]$ifpid, sum)/tapply(tmpwt, newdat$ifpdat[[j]]$ifpid, sum)
                mmde <- sum(tmpwt * rowbriers)
                decomp <- c(summary(mde), sd(mde)) #c(2*mmde, mmde, 2*mmde, 2*mmde, mmde, 2*mmde)
                reslist[[j]][i,] <- decomp
            }
        }
    }
    components <- sapply(reslist, apply, 2, mean)

    if(scale){
        ## run once more to get mmde
        sysmmde <- rep(NA, ngrp)
        for(j in 1:ngrp){
            tmpf <- shuffdat$fdat[[j]]
            tmpd <- shuffdat$ddat[[j]]
            tmpbin <- shuffdat$newbins[[j]]
            nbin <- shuffdat$nbin
            tmpwt <- newdat$wtdat[[j]]

            decomp <- calcDecomp(tmpf, tmpd, tmpbin, tmpwt, nbin, scale = FALSE)
            decomp <- decomp$res

            sysmmde[j] <- -decomp[1] + decomp[2] + decomp[6]
        }
        components <- rbind(components, sysmmde)
    } else {
        components <- rbind(components, apply(components, 2, function(x){
          -x[1] + x[2] + x[6]}))
    }

    colnames(components) <- grpnames
    if(bin){
        rownames(components) <- c("discrim", "miscal", "deltaf",
                                  "miscal_lg", "cov", "unc", "brier")
    }
    if(scale) rownames(components)[c(2,4)] <- c("cal", "cal_lg")
  
    if(bin){
        components <- rbind(components,
                            as.matrix(origmmde[c(4,7,3,1:2,5:6),],
                                      nrow(origmmde), ncol(origmmde)))
        rownames(components)[8:14] <- c("orig_brier","mde_sd","mde_median","mde_min","mde_q25","mde_q75","mde_max")
        print(round(components, 3))
        names(reslist) <- names(qcalib) <- names(qdiscrim) <- grpnames
        res <- list(components = components, resamples = reslist,
                    qcalib = qcalib, qdiscrim = qdiscrim, scale = scale)
    } else {
        res <- components
    }

    res
}



## test code
if(FALSE){
  f <- round(matrix(runif(20, 0, .2), 10, 2), 1)
  f <- cbind(f, rep(NA, 10))
  f[,3] <- 1 - apply(f,1,sum,na.rm=T)
  bin <- as.numeric(factor(apply(f,1,paste,collapse="")))
  wt <- rep(.1,10)

  d <- t(rmultinom(10, 1, c(.3,.3,.4)))
  dvec <- apply(d, 1, function(x) which(x==1))

  tmp <- bdecomp(f, dvec, wt=wt)


  ## old code:
  source("~/tmp/mmdeDecomp0.5/mmdeDecomp.R")
  source("~/tmp/mmdeDecomp0.5/expandOrd.R")
  source("~/tmp/mmdeDecomp0.5/calcDecomp.R")
  source("~/tmp/mmdeDecomp0.5/shuffle.R")
  source("~/tmp/mmdeDecomp0.5/mmdeIntervals.R")

  set.seed(100)
  res <- mmdeDecomp("~/tmp/mmdeDecomp0.5/acedat.csv", "~/tmp/mmdeDecomp0.5/qtype.csv", niter=50, debug=TRUE, scale=FALSE)
  mmdeIntervals(res)

  ## new code
  load_all('~/Projects/scoring/pkg')
  fs <- read.csv("~/tmp/mmdeDecomp0.5/acedat.csv")
  qs <- read.csv("~/tmp/mmdeDecomp0.5/qtype.csv")

  nday <- with(subset(fs, system==1), tapply(outcome, ifpid, length))
  nq <- length(unique(fs$ifpid))
  fs$nday <- nday[match(fs$ifpid, as.numeric(names(nday)))]
  fs$nq <- nq
  fs$wt <- 1/(fs$nday * fs$nq)

  ## matches fairly well
  set.seed(100)
  newres <- bdecomp(fs[,1:5], fs$outcome, fs$system, wt=fs$wt,
                    qid=fs$ifpid, qtype=qs, resamples=50)
  mmdeIntervals(newres)
}
