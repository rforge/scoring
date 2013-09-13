calcscore <- function(object, ...){
    UseMethod("calcscore")
}

calcscore.formula <- function(object, fam="pow", param, data, bounds=NULL, reverse=FALSE, ...){

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
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    outcome <- model.response(mf, "any")
    forecast <- mf[,2:ncol(mf)]
    if(missing(param)) param <- c(2, rep(1/max(2,NCOL(forecast)), max(2,NCOL(forecast))))

    do.call("calcscore.default", list(object=forecast, outcome=outcome, fam=fam, param=param, bounds=bounds, reverse=reverse))

}
            
calcscore.default <-
function(object, outcome, fam="pow", param=c(2,rep(1/max(2,NCOL(forecast)),max(2,NCOL(forecast)))), bounds=NULL, reverse=FALSE, ...){

    ## For deprecated scaling argument
    dots <- list(...)

    if("scaling" %in% names(dots)){
      if(dots$scaling) bounds <- c(0,1)
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
      ss <- rowSums(forecast)
      if(any(ss != 1)){
        warning("Forecasts in some rows do not sum to 1; they were scaled to sum to 1.")
      }
    }
    ## Check beta family params
    if(fam=="beta"){
        if(any(param <= -1)) stop("Beta family parameters must be greater than -1")
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
    npars <- length(param)
    
    ## For fam=pow or sph, check to ensure that baseline params sum to 1.
    if(fam %in% c("pow","sph")){
        if(sum(param[2:npars]) != 1){
            ## If two alternatives and one baseline, take
            ## complement
            if(npars==2){
                param <- c(param, 1-param[2])
            } else {
                if(npars != (nalts+1)) stop("Length of param is incorrect.\n")
                warning("Baseline parameters were scaled to sum to 1.")
                param[2:npars] <- param[2:npars]/sum(param[2:npars])
            }
        }
    }    
    ## END ERROR CHECKING

    ## Create data matrix
    if(NCOL(forecast)==1) forecast <- cbind(1-forecast, forecast)
    datmat <- cbind(forecast, outcome)
    ## Obtain unscaled scores
    sc <- scoreitems(param, datmat, fam)

    ## Scale if desired
    if(!is.null(bounds) | reverse){
        ## Aug 26 2013: scoreitems appears to handle multiple alts
        ## 2-alternative examples yield same results as before
        scalefactor <- scalescores(param, fam)

        lbound <- ifelse(is.na(bounds[1]), 0, bounds[1])
        ubound <- ifelse(is.na(bounds[2]), 1, bounds[2])

        ## Note: reverse is logical but used as 0/1 below.
        if(fam=="beta"){
            ## TODO consider other scaling for beta family?
            sc <- sc/scalefactor
        } else {
            sc <- (sc - scalefactor[1])/diff(scalefactor)
        }
        if(reverse) sc <- 1 - sc
        sc <- lbound + (ubound - lbound)*sc
    }
    if(any(is.na(sc))){
        stop("Problem with score calculation.  Ensure parameter values are valid.")
    }
    sc
}
