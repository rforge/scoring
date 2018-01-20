scoreitems <- function(param, data, fam, ordered, decomp=FALSE, group=NULL, decargs=NULL){
    ordmix <- (length(ordered) > 1) & (!all(ordered == ordered[1]))
    if(ordmix){
        ordrows <- which(ordered)
        urows <- which(!ordered)
    } else if(length(ordered) > 0 & ordered[1]){
        ordrows <- 1:nrow(data)
        urows <- NULL
    } else {
        ordrows <- NULL
        urows <- 1:nrow(data)
    }

    item.res <- rep(NA, nrow(data))
    if(length(ordrows) > 0){
        ## deal with differing numbers of alternatives
        nalts <- apply(!is.na(data[ordrows,1:(ncol(data)-1)]), 1,
                       sum)
        totalts <- unique(nalts)

        for(i in 1:length(totalts)){
            tmprow <- which(nalts == totalts[i])

            parm <- param
            if(nrow(param) > 1) parm <- as.matrix(param[ordrows[tmprow],])
            tmpfam <- fam
            if(length(fam) > 1) tmpfam <- fam[ordrows[tmprow]]
            item.res[ordrows[tmprow]] <- ordwrap(data[ordrows[tmprow],c(1:totalts[i],ncol(data))],
                                                 parm,
                                                 tmpfam)
        }
    }
    if(length(urows) > 0){
        nalts <- ncol(data)-1

        parm <- param
        if(nrow(param) > 1) parm <- as.matrix(param[urows,])
        tmpfam <- fam
        if(length(fam) > 1) tmpfam <- fam[urows]
        scorefun <- paste(tmpfam,"fam",sep="")

        datau <- cbind.data.frame(data[urows,], scorefun, parm)

        item.res[urows] <- apply(datau, 1, function(x){
            scfun <- x[(nalts+2)]
            ## allow differing numbers of alts per row:
            fvals <- which(!is.na(x[1:nalts]))
            outs <- match(x[(nalts+1)], fvals)
            x <- as.numeric(x[-(nalts+2)])
            do.call(scfun, list(x[fvals], x[(nalts+1)],
                                x[(nalts+2):length(x)], scfun))
        })
    }

    if(decomp){
        decargs <- c(decargs, list(forecast=data[,1:nalts],
                                   outcome=data[,(nalts+1)],
                                   group=group))
        d.res <- do.call("bdecomp", decargs)

        item.res <- list(rawscores=item.res, decomp=d.res)
    }
    
    item.res
}
