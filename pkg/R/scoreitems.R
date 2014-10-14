scoreitems <- function(param, data, fam, ordered, decomp=FALSE, group=NULL){
    if(ordered){
        item.res <- ordwrap(data, param, fam)
    } else {
        scorefun <- paste(fam,"fam",sep="")

        nalts <- ncol(data)-1

        data <- cbind.data.frame(data, scorefun, param)

        item.res <- apply(data, 1, function(x){
            scfun <- x[(nalts+2)]
            x <- as.numeric(x[-(nalts+2)])
            do.call(scfun, list(x[1:nalts], x[(nalts+1)],
                                x[(nalts+2):length(x)], scfun))
        })
    }

    if(decomp){
        if(nalts>2) stop("Brier score decompositions are only valid for 2-alternative forecasts.\n")
        d.res <- bdecomp(data[,1:nalts], data[,(nalts+1)], group)

        item.res <- list(rawscores=item.res, decomp=d.res)
    }
    
    item.res
}
