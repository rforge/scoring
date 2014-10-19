"glmscore" <- function(y, X, fam, param)
    {
        discrep <- function(beta, y, X, fam, param){
            preds <- plogis(X %*% beta)

            res <- calcscore(y ~ preds, fam = fam, param = param,
                             bounds = c(0, 1))

            mean(res)
        }

        betastart <- c(qlogis(mean(y)), rep(0, (ncol(X)-1)))

        if(fam=="beta"){
            grad <- betagrad
        } else {
            grad <- NULL
        }
        beta <- nlminb(betastart, discrep, y=y, X=X, fam=fam,
                       param=param, lower=-10, upper=10,
                       gradient=grad)

        beta
    }

betagrad <- function(beta, y, X, param, ...){
    preds <- plogis(X %*% beta)
    tmp <- (y - preds) * preds^(param[1]-1) *
           (1 - preds)^(param[2]-1) *
           binomial()$mu.eta(X %*% beta)

    -cbind(mean(tmp * X[,1]), mean(tmp * X[,2]))
}


if(FALSE){
    library(scoring)
    library(smdata)
    data(finance)
    source('glmscore.R')

    m1 <- glm(corr ~ probc, data=finance, family=binomial)

    ## warnings from bounds argument
    m1b <- with(finance, glmscore(corr,
                                  cbind(rep(1, length(corr)),
                                        probc),
                                  fam = "beta", param=c(0,0)))

    ## now brier score
    m1c <- with(finance, glmscore(corr,
                                  cbind(rep(1, length(corr)),
                                        probc),
                                  fam = "beta", param=c(1,1)))

    m1d <- with(finance, glmscore(corr,
                                  cbind(rep(1, length(corr)),
                                        probc),
                                  fam = "beta", param=c(11,15)))

    ## Other family
    m1e <- with(finance, glmscore(corr,
                                  cbind(rep(1, length(corr)),
                                        probc),
                                  fam = "pow", param=c(4, .1)))
    
    ## family could work like this, but doesn't seem
    ## tenable because this is not an exponential family
    ## distribution
    scfam <- binomial()
    scfam$family <- "scoring"
    scfam$dev.resids <- function(y, mu, wt){
        if(length(mu)==1) mu <- rep(mu, length(y))
        wt * calcscore(y ~ mu, fam="pow", param=c(0, 0), bounds=c(0,1))
    }
    m1f <- glm(corr ~ probc, data=finance, family=scfam)
}
