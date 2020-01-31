#' sch
#'
#' build the knn indices based on schlittgens' algorithm
#'
#'
#' @param X the design matrix
#' @param y the response vector
#' @param rs number of repetitions
#' @param k number of clusters
#' @return Returns an n by k matrix; each row a sorted data indices
#' @export sch

sch <- function(X,y,rs,k){

    n <- nrow(X)
    p <- ncol(X)
    nbr <- matrix(0, nrow = n, ncol = n)
    member.matrix <- matrix(0, nrow = rs, ncol = n)

    for (r in c(1:rs)){


        w1 <- c((1-(k-1)*(k*2)^-1),((k*2)^-1*rep(1,k-1)))

        w <- matrix(0,n,k)

        for (i in c(1:n)){
            w[i,] <- sample(w1,k,replace=FALSE) }

        beta <- matrix(0,p,k)
        res <- matrix(0,n,k)

        it <- 0
        betaold <- matrix(1,p,k)


        while ((it<100) & ( sum(abs(beta-betaold))> 0.00001 ) ){

            it <- it+1
            betaold <- beta

            for (j in c(1:k)){
                out <- lm.wfit(X, y, w[,j])
                res[,j] <- out$residuals
                beta[,j] <- out$coefficients
            }

            ares <- res^-2 #square of residual inverse
            w <- ares/apply(ares,1,sum)

        }

        member <- apply(w,1,max)
        index <- t(matrix(c(1:k),k,n))
        member <- apply(index*(w==member),1,sum)


        for (i in 1:n){
            for (j in 1:n){
                if(member[i] == member[j]){
                    nbr[i,j] <- nbr[i,j] + 1
                }
            }
        }
    }

    nbr.index <- t(sortIndex(nbr, decrease = T))
    nbr.index1 <- matrix(0, nrow = n, ncol = n-1)

    for (i in 1:n){
        me.index <- which(nbr.index[i,] == i)
        nbr.index1[i,] <- nbr.index[i,-me.index]

    }

    nbr.index1
}


#' clusterWise
#'
#' build the knn indices based on schlittgens' algorithm
#'
#' @param x the design matrix (intercept term should be included if any)
#' @param y the response vector
#' @param k number of clusters
#' @param rep number of replications; default  = 30
#' @return Returns the membership labels
#' @export clusterWise

clusterWise <- function(x,y,k,rep = 30){

    n <- nrow(x)
    p <- ncol(x)
    ind <- c(1:n)
    w1 <- c((1-(k-1)*(k*2)^-1),((k*2)^-1*rep(1,k-1)))

    beta.rep <- array(1,dim=c(p,k,rep))
    member.rep <-  matrix(0,n,rep)
    R2.rep <-  rep(NA,rep)
    AIC.rep <-  rep(NA,rep)

    for (W in c(1:rep)){

        w <- matrix(0,n,k)

        for (i in c(1:n)){
            w[i,] <- sample(w1,k,replace=FALSE) }

        beta <- matrix(0,p,k)
        res <- matrix(0,n,k)
        it <- 0
        betaold <- matrix(1,p,k)

        while ((it<100) & ( sum(abs(beta-betaold))> 0.00001 ) ){
            it <- it+1
            betaold <- beta

            for (j in c(1:k)){
                out <- lm.wfit(x, y, w[,j])
                res[,j] <- out$residuals
                beta[,j] <- out$coefficients
            }

            ares <- res^-2
            w <- ares/apply(ares,1,sum)
        }

        member <- apply(w,1,max)
        ed <- t(matrix(c(1:k),k,n))
        member <- apply(ed*(w==member),1,sum)

        w <- rep(0,n)
        AIC <- 0
        zaehl <- rep(0,k)

        for (j in c(1:k)){

            indj <- ind[member==j]
            nj <- length(indj)

            if (nj > p+1){
                out <- lm.fit(x[indj,], y[indj])
                beta[,j] <- out$coefficients
                s <- var(out$residuals)*(nj-1)/(nj-p)
                w[indj] <- 1/s
                zaehl[j] <-  sum(y[indj]*(x[indj,]%*%beta[,j]))/s
                AIC <- AIC + log(s)*nj
            }

            else{beta[,j] <- rep(0,p)}

        }

        R2 <- (sum(zaehl) - (sum(y*w)^2)/sum(w))/( sum((y^2)*w) - (sum(y*w)^2)/sum(w))

        beta.rep[,,W] <- beta
        member.rep[,W] <- member
        R2.rep[W] <-  R2
        AIC.rep[W] <-  AIC
    }

    o <- order(R2.rep)
    R2 <- R2.rep[o[rep]]
    member <- member.rep[,o[rep]]
    beta <- beta.rep[,,o[rep]]
    beta <- as.matrix(beta)
    R2adj <- 1-(1-R2)*(n-1)/(n-k*p)
    AIC <- AIC.rep[o[rep]] + 2*k*p
    BIC <- AIC.rep[o[rep]] + log(n)*k*p
    member
}




