"glmscore" <- function(y, X, fam, param, usehess=FALSE)
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
            if(usehess){
              hess <- betahess
            } else {
              hess <- NULL
            }
        } else {
            grad <- NULL
        }
        beta <- nlminb(betastart, discrep, y=y, X=X, fam=fam,
                       param=param, gradient=grad, hessian=hess)

        beta
    }

betagrad <- function(beta, y, X, fam, param){
    preds <- plogis(X %*% beta)
    tmp <- (y - preds) * preds^(param[1]-1) *
           (1 - preds)^(param[2]-1) *
           binomial()$mu.eta(X %*% beta)

    -cbind(mean(tmp * X[,1]), mean(tmp * X[,2]))
}

betahess <- function(beta, y, X, fam, param){
    preds <- plogis(X %*% beta)
    omega <- preds^(param[1] - 1) * (1 - preds)^(param[2] -1)

    hess <- matrix(0, ncol(X), ncol(X))
    tmp <- omega * binomial()$mu.eta(X %*% beta)^2 -
           (y - preds) * (param[1]*(1 - preds) - param[2]*preds)

    for(i in 1:nrow(X)){
      hess <- hess + tmp[i] * X[i,] %*% t(X[i,])
    }

    hess/nrow(X)
}

## not yet working for non-1 baseline parameters
powgrad <- function(beta, y, X, fam, param){
    ## first entry of q is for d=0, second entry for d=1
    if(length(param)==1){
      q <- c(1, 1)
    } else {
      q <- param[2:3]
    }

    preds <- plogis(X %*% beta)
    pred1 <- (1 - y) - (1 - 2*y)*preds
    dlogi <- binomial()$mu.eta(X %*% beta)

    tmp1 <- -(1 - 2*y) * ((pred1)/q[(y+1)])^(param[1]-2) * (dlogi/q[(y+1)])
    tmp2 <- ((preds)/q[2])^(param[1]-1) * (dlogi/q[2])
    tmp3 <- ((1-preds)/q[1])^(param[1]-1) * (dlogi/q[1])

    tmp4 <- tmp1 - (tmp2 - tmp3)
    grad <- t(tmp4) %*% X

    grad/nrow(X)
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

    X <- cbind(rep(1, length(finance$corr)), finance$probc)
    m1bh <- betahess(m1b$par, finance$corr, X, c(0,0))
    m1bvc <- solve(nrow(X) * m1bh)
    m1vc <- vcov(m1)
    all.equal(m1bvc, m1vc)
    
    ## now brier score
    m1c <- with(finance, glmscore(corr, X, fam = "beta", param=c(1,1),
                                  usehess=FALSE))
    ## S.E. for intercept blows up
    solve(nrow(X) * betahess(m1c$par, finance$corr, X, NULL, c(1,1)))

    ## match from power family?
    m1c2 <- with(finance, glmscore(corr, X, fam = "pow", param=2))

    ## strange score, difficult to converge due to flat parameter space
    ## (and usehess=FALSE falsely converges)
    m1d <- with(finance, glmscore(corr, X, fam = "beta", param=c(11,15), usehess=TRUE))

    m1d <- with(finance, glmscore(corr, X, fam="beta", param=c(2.25, 3.75), usehess=TRUE))
    
    ## log score from power family
    m1e <- with(finance, glmscore(corr, X, fam="pow", param=c(2, .2, .8)))

    
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
