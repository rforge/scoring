scalescores <-
function(pars, fam="pow"){
  ps <- c(0,1)

  ## Shortcut for scaling scores for >2 alts
  ## Min must be achieved when forecast 1 for lowest baseline
  ##     and lowest baseline occurs.
  ## Max must be achieved when forecast 0 for highest baseline
  ##     and 1 for lowest and highest baseline occurs.
  
  if(fam=="beta"){
      xplier <- max(betafam(c(1,0), d=2, param=pars),
                    betafam(c(0,1), d=1, param=pars))

      if(xplier == Inf){
          warning("Scaling does not work because maximum possible score is Inf.")
          xplier <- 1
      }
      
  } else if(fam=="pow" | fam=="sph"){
      nalts <- length(pars)-1
      ## Baseline parameters
      baselines <- pars[2:(nalts+1)]
      ## Which are largest and smallest?
      maxbase <- which(baselines==max(baselines))[1]
      minbase <- which(baselines==min(baselines))
      ## Take last category, in case they are all equal
      minbase <- minbase[length(minbase)]
      
      fore <- out <- rep(0,nalts)
      fore[minbase] <- 1
      out[maxbase] <- 1
      tmpsc <- calcscore(c(minbase,minbase) ~ rbind(fore,out), fam=fam, param=pars)
      ##scmin <- calcscore(matrix(fore,1,nalts) ~ matrix(fore,1,nalts), fam=fam,
      ##                   param=pars)
      ##scmax <- calcscore(matrix(fore,1,nalts) ~ matrix(out,1,nalts), fam=fam,
      ##                   param=pars)
      
      xplier <- tmpsc #c(scmin, scmax)
      if(xplier[2]==Inf){
          warning("Scaling does not work because maximum possible score is Inf.")
          xplier <- c(0,1)
      }
  }

  xplier
}

if(FALSE){
    ## Proof that above yields min and max values for 3 alts:
    p <- seq(0,1,.01)
    y <- expand.grid(p,p,p)
    ysum <- apply(y,1,sum)
    y <- y[ysum==1,]
    out1 <- t(matrix(c(1,0,0),3,nrow(y)))
    out2 <- t(matrix(c(0,1,0),3,nrow(y)))
    out3 <- t(matrix(c(0,0,1),3,nrow(y)))

    sc1 <- calcscore(y ~ out1,fam="pow",param=c(3,.1,.5,.4))
    sc2 <- calcscore(y ~ out2,fam="pow",param=c(3,.1,.5,.4))
    sc3 <- calcscore(y ~ out3,fam="pow",param=c(3,.1,.5,.4))

    c(max(sc1),max(sc2),max(sc3))
    y[which(sc2==max(sc2)),]
    c(min(sc1),min(sc2),min(sc3))
    y[which(sc1==min(sc1)),]
}
