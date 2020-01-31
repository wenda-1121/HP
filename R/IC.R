#' log.factorial
#'
#' calculate c( log{(n_1)!}, log{(n_2)!}, ... log{(n_k)!} )
#'
#' @param n.k a vector: c(n_1, n_2, ... n_k) where the log-factorial is applied to.
#' @return c( log{(n_1)!}, log{(n_2)!}, ... log{(n_k)!} )
#' @noRd

log.factorial <- function(n.k){

    res <- rep(0,length(n.k))

    for (i in 1:length(n.k)){
        s <- 1:n.k[i]
        res[i] <- sum(log(s))
    }

    return(res)

}


#' ic
#'
#' calculate the aic/bic/mdl of a given clusterwise regression model
#'
#' @param X the design matrix
#' @param y the response vector
#' @param membership a vector that records the membership of each data case
#' @param type "aic", "bic" or "mdl"
#' @return Returns the specified information criterion
#' @export ic

ic <- function(X,y, membership, type = "mdl"){

    p <- ncol(X)
    n <- nrow(X)

    K <- length(unique(membership))
    n.k <- rep(0, K)

    for (i in 1:K){

        index.Gk <- which(membership == i)
        size.k <- length(index.Gk)

        if (size.k <= p){

            return(NaN)
            print("sample size too small for some group(s)!")

        }

        n.k[i] <- size.k

    }

    sigma2.hat <- rep(0, K)

    for (i in 1:K){

        index.Gk <- which(membership == i)

        y.Gk <- y[index.Gk]
        X.Gk <- X[index.Gk,]
        sigma2.hat[i] <- mean( (residuals(lm(y.Gk ~ X.Gk - 1)))^2 )
    }

    alpha <- n.k / n
    term1 <- sum(n.k * log(sigma2.hat))


    if (type == "bic"){

        ic <- term1 + K * (p + 1) * log(n)
    }

    if (type == "aic"){

        ic <- term1 +  2 * K * (p + 1)
    }


    if (type == "mdl"){

        ic <- term1 +  (p+1) * sum(log(n.k)) +
            2 * sum(log.factorial(n.k)) -
            2 * sum(n.k * log(alpha)) +  (K-1) * log(n)

    }

    return(ic)
}
