scoreitems <- function(param, data, fam, ordered){
    if(ordered){
        item.res <- ordwrap(data, param, fam)
    } else {
        scorefun <- paste(fam,"fam",sep="")

        nalts <- ncol(data)-1

        item.res <- apply(data, 1, function(x){
            do.call(scorefun, list(x[1:nalts], x[length(x)], param, fam))
        })
    }
    
    item.res
}
