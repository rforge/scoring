"glmscore" <- function(y, X, fam, param, usehess=FALSE)
    {
        discrep <- function(beta, y, X, fam, param){
            preds <- plogis(X %*% beta)

            res <- calcscore(y ~ preds, fam = fam, param = param, bounds=c(0,1e3))

            mean(res)
        }

        betastart <- coef(glm.fit(X, y, family=binomial()))

        if(fam=="beta"){
            grad <- betagrad
            if(usehess){
              hess <- betahess
            } else {
              hess <- NULL
            }
        } else if(fam=="pow"){
            grad <- powgrad
            if(usehess){
              hess <- powhess
            } else {
              hess <- NULL
            }            
        }
        ## was nlminb
        beta <- optim(betastart, discrep, gr=grad, y=y, X=X, fam=fam,
                       param=param, hessian=TRUE, method="BFGS") #gradient=grad, hessian=hess)

        beta
    }

betagrad <- function(beta, y, X, fam, param, sc=FALSE){
    preds <- plogis(X %*% beta)
    tmp <- (y - preds) * preds^(param[1]-1) *
           (1 - preds)^(param[2]-1) *
           binomial()$mu.eta(X %*% beta)

    if(sc) return(-cbind(tmp * X[,1], tmp * X[,2]))
    -cbind(mean(tmp * X[,1]), mean(tmp * X[,2]))
}

betahess <- function(beta, y, X, fam, param){
    preds <- plogis(X %*% beta)
    omega <- preds^(param[1] - 1) * (1 - preds)^(param[2] - 1)

    hess <- matrix(0, ncol(X), ncol(X))
    tmp <- omega * binomial()$mu.eta(X %*% beta)^2 -
           (y - preds) * (param[1] * preds^(param[1]-1) * (1-preds)^(param[2]) +
                          param[2] * preds^(param[1]) * (1-preds)^(param[2]-1))

    for(i in 1:nrow(X)){
      hess <- hess + tmp[i] * tcrossprod(X[i,])
    }

    hess/nrow(X)
}

betapsi <- function(beta, y, X, fam, param){
    preds <- plogis(X %*% beta)
    psi <- -(y - preds) * preds^(param[1] - 1) * (1 - preds)^(param[2] - 1)
    psiprime <- (param[1] - y) * preds^(param[1] - 1 - y) *
                (1 - preds)^(param[2] - 1 + y) +
                (param[2] - 1 + y) * preds^(param[1] - y) * (1 - preds)^(param[2] - 2 + y)
    list(psi = psi, psiprime = psiprime)
  }

betavcov <- function(beta, y, X, fam, param){
    tmp <- betagrad(beta, y, X, fam, param, sc=TRUE)
    tmp2 <- betahess(beta, y, X, fam, param)

    (1/nrow(X)) * solve(tmp2) %*% cov(tmp) %*% t(solve(tmp2))
  }

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

    tmp1 <- -(1 - 2*y) * q[(y+1)]^(-(param[1]-1)) * pred1^(param[1]-2)
    tmp2 <- ((preds)/q[2])^(param[1]-1)
    tmp3 <- ((1-preds)/q[1])^(param[1]-1)

    tmp4 <- (tmp1 - (tmp2 - tmp3)) * dlogi
    grad <- t(tmp4) %*% X

    ## negative because this is how family is defined in scoring package
    -grad/nrow(X)
}

powhess <- function(beta, y, X, fam, param){
    if(length(param)==1){
      q <- c(1, 1)
    } else {
      q <- param[2:3]
    }

    preds <- plogis(X %*% beta)
    pred1 <- (1 - y) - (1 - 2*y)*preds
    dlogi <- binomial()$mu.eta(X %*% beta)

    tmp1 <- (1 - 2*y)^2 * (param[1]-2) * pred1^(param[1]-3) * q[(y+1)]^(-(param[1]-1))
    tmp2 <- (param[1]-1) * preds^(param[1]-2) * q[2]^(-(param[1]-1))
    tmp3 <- (param[1]-1) * (1 - preds)^(param[1]-2) * q[1]^(-(param[1]-1))

    tmp4 <- (tmp1 - (tmp2 + tmp3)) * dlogi^2

    hess <- matrix(0, ncol(X), ncol(X))

    for(i in 1:nrow(X)){
      hess <- hess + tmp4[i] * tcrossprod(X[i,])
    }

    ## negative because this is how family is defined in scoring package
    -hess/nrow(X)
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
    m1bh <- betahess(m1b$par, finance$corr, X, "beta", c(0,0))
    m1bvc <- solve(nrow(X) * m1bh)
    m1vc <- vcov(m1)
    all.equal(m1bvc, m1vc)
    tmp <- betavcov(m1b$par, finance$corr, X, "beta", c(0,0))
    
    ## now brier score
    m1c <- with(finance, glmscore(corr, X, fam = "beta", param=c(1,1),
                                  usehess=FALSE))
    ## S.E. for intercept blows up this way:
    solve(nrow(X) * betahess(m1c$par, finance$corr, X, NULL, c(1,1)))
    ## but that's because it was wrong
    betavcov(m1c$par, finance$corr, X, "beta", c(1,1))

    ## match from power family?
    m1c2 <- with(finance, glmscore(corr, X, fam = "pow", param=2))
    powgrad(m1c2$par, finance$corr, X, "pow", 2)
    ## not sure why this is smaller than the one from m1c
    solve(nrow(X) * powhess(m1c2$par, finance$corr, X, NULL, 2))
    
    ## strange score, difficult to converge due to flat parameter space
    ## (and usehess=FALSE falsely converges)
    m1d <- with(finance, glmscore(corr, X, fam = "beta", param=c(11,15), usehess=TRUE))
    betavcov(m1d$par, finance$corr, X, "beta", c(11,15))

    m1d <- with(finance, glmscore(corr, X, fam="beta", param=c(2.25, 3.75), usehess=TRUE))
    betavcov(m1d$par, finance$corr, X, "beta", c(2.25,3.75))
    
    ## log score from power family
    m1e <- with(finance, glmscore(corr, X, fam="pow", param=c(4.3, .1, .9)))
    powgrad(m1e$par, finance$corr, X, "pow", c(4.3, .1, .9))
    powhess(m1e$par, finance$corr, X, "pow", c(4.3, .1, .9))
}
