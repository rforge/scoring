calcscore <- function(object, ...){
    UseMethod("calcscore")
}

calcscore.formula <- function(object, fam="pow", param, data, bounds=NULL, reverse=FALSE, ordered=FALSE, ...){

    ## For deprecated scaling argument
    dots <- list(...)

    if("scaling" %in% names(dots)){
      if(dots$scaling) bounds <- c(0,1)
    }

    if(missing(data)) data <- environment(object)
    mf <- match.call()
    m <- match(c("object", "data"), names(mf), 0L)
    names(mf)[m[1L]] <- "formula"
    mf <- mf[c(1L, m)]

    ## Get outcomes and forecasts
    mf[[1L]] <- as.name("model.frame")
    ## to send rows with NA through:
    if(length(mf) == 3L){
      mf[[4L]] <- na.pass
      names(mf)[4] <- "na.action"
      #attr(mf[[3L]], 'na.action') <- na.pass
    }
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    outcome <- model.response(mf, "any")
    forecast <- mf[,2:ncol(mf)]
    if(missing(param)) param <- c(2, rep(1/max(2,NCOL(forecast)), max(2,NCOL(forecast))))

    do.call("calcscore.default", list(object=forecast, outcome=outcome, fam=fam, param=param, bounds=bounds, reverse=reverse, ordered=ordered, ...=...))

}
            
calcscore.default <-
function(object, outcome, fam="pow", param=c(2,rep(1/max(2,NCOL(forecast)),max(2,NCOL(forecast)))), bounds=NULL, reverse=FALSE, ordered=FALSE, ...){

    dots <- list(...)

    ## For deprecated scaling argument
    if("scaling" %in% names(dots)){
      if(dots$scaling) bounds <- c(0,1)
    }

    ## For Brier score decompositions
    if("decomp" %in% names(dots)){
        decomp <- dots$decomp
        if("group" %in% names(dots)){
            group <- dots$group
        } else {
            group <- rep(1, NROW(object))
        }

        if("decompControl" %in% names(dots)){
            decargs <- dots$decompControl
        } else {
            decargs <- NULL
        }
    } else {
        decomp <- FALSE
        group <- 1
    }

    ## Error checks:
    ## Make sure outcome is a vector, then convert it to numeric
    if(!(NROW(outcome)==1 | NCOL(outcome)==1)) stop("outcome must be a vector, not a matrix.\n")

    forecast <- object
    outcome <- as.numeric(outcome)
    noutcomes <- length(outcome)
    ## Make sure number of forecasts == number of outcomes
    if(NROW(forecast) != noutcomes) stop("Number of forecasts does not match number of outcomes. Check dimension of forecast.\n")

    ## If it is two alternatives and forecasts for only one outcome
    ## are supplied, convert to numeric
    if(NCOL(forecast) == 1) forecast <- as.numeric(forecast)
    
    ## Make sure beta family is not used with >2 alternatives.
    if(NCOL(forecast)>2 & fam=="beta") stop("Beta family cannot be used with >2 alternatives.\n")
    ## Check that forecast rows sum to 1.  If not, rescale and warn.
    if(is.matrix(forecast) | is.data.frame(forecast)){
      ss <- rowSums(forecast, na.rm=TRUE)
      if(any(ss != 1)){
        warning("Forecasts in some rows do not sum to 1; they were scaled to sum to 1.")
      }
    }
    ## Check beta family params
    if(fam=="beta"){
        if(any(param <= -1)) stop("Beta family parameters must be greater than -1")
        if(any(ordered)){
            ordered <- FALSE
            warning("ordered=TRUE has no impact on beta family scores.")
        }
    }
    ## If forecast is numeric, add 1 to outcome
    if(NCOL(forecast)==1){
        ## This is a binary problem; assume outcomes are 0/1
        ## and add 1
        ## (place outcome==1 forecasts in second column)
        outcome <- outcome+1
    }
    ## Assume outcome is vector, so don't need:
    ## if(NCOL(outcome) > 1) outcome <- apply(outcome, 1, function(x) which(x==1))

    ## Number of alternatives and number of parameters
    nalts <- ifelse(NCOL(forecast) <= 2, 2, NCOL(forecast))
    if(class(param)=="numeric") param <- matrix(param, 1, length(param))
    npars <- NCOL(param)
    nparrows <- NROW(param)

    ## For fam=pow or sph, check to ensure that baseline params sum to 1.
    if(fam %in% c("pow","sph")){
        ## If length(param)==1, then assume this is a rule without
        ## a baseline.
        if(npars > 1 & nparrows == 1){
            if(sum(param[1, 2:npars]) != 1){
                ## If two alternatives and one baseline, take
                ## complement
                if(npars==2){
                    if(NCOL(forecast)==1){
                        param <- matrix(c(param[1,1], 1-param[1,2], param[1,2]),
                                        1, 3)
                    } else {
                        warning("Only one baseline parameter supplied. This parameter is assumed to correspond to the alternative associated with the first column of forecasts.\n")
                        param <- matrix(c(param, 1-param[1,2]), 1, 3)
                    }
                } else {
                    if(npars != (nalts+1)) stop("Length of param is incorrect.\n")
                    warning("Baseline parameters were scaled to sum to 1.")
                    param[1, 2:npars] <- param[1, 2:npars]/sum(param[1, 2:npars])
                }  ## npars==2
            }  ## sum != 1
            param <- t(matrix(param, ncol(param), noutcomes))
            fam <- rep(fam, noutcomes)
        } else if(npars == 1 & nparrows == 1){
            param <- matrix(rep(param, noutcomes), noutcomes, 1)
            fam <- rep(fam, noutcomes)
        } else if(npars > 1 & nparrows == noutcomes){
            ## Ensure it is a matrix
            param <- as.matrix(param)
            if(length(fam)==1) fam <- rep(fam, noutcomes)
        } else {
            stop("Problem with param matrix.  Check that number of param rows equal number of forecast rows.")
        }
    }
    ## END ERROR CHECKING

    ## Create data matrix
    if(NCOL(forecast)==1) forecast <- cbind(1-forecast, forecast)
    datmat <- cbind(forecast, outcome)
    ## Obtain unscaled scores
    sc <- scoreitems(param, datmat, fam, ordered, decomp, group,
                     decargs)

    ## If decomp=TRUE, sc is a list:
    if(decomp){
        sclist <- sc
        sc <- sc$rawscores
    }

    ## Scale if desired
    if(!is.null(bounds)){
        ## TODO Apply scalescores() to each unique scoring rule
        if(nparrows > 1 | !all(ordered == ordered[1])){
            warning("Scaling currently unavailable when the scoring rule differs by row.")
        } else {
            scalefactor <- scalescores(param[1,], fam[1], ordered, nalts)

            lbound <- ifelse(is.na(bounds[1]), 0, bounds[1])
            ubound <- ifelse(is.na(bounds[2]), 1 + lbound, bounds[2])

            if(fam[1]=="beta"){
                sc <- sc/scalefactor
            } else {
                sc <- (sc - scalefactor[1])/diff(scalefactor)
            }
            sc <- lbound + (ubound - lbound)*sc

            if(reverse){
                sc <- lbound + ubound - sc
            }
        }
    }

    if(reverse & is.null(bounds)) sc <- -sc

    if(any(is.na(sc))){
        warning("Some scores are NA. This may be due to missing data in your forecasts or outcomes, or an ill-defined param argument.")
    }

    if(decomp){
        sclist$rawscores <- sc
        sc <- sclist
    }
    
    sc
}
