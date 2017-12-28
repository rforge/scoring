scoreitems <- function(param, data, fam, ordered, decomp=FALSE, group=NULL, decargs=NULL){
    ## FIXME allow ordered to vary by row
    if(ordered){
        item.res <- ordwrap(data, param, fam)
    } else {
        scorefun <- paste(fam,"fam",sep="")

        nalts <- ncol(data)-1

        data <- cbind.data.frame(data, scorefun, param)

        item.res <- apply(data, 1, function(x){
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
