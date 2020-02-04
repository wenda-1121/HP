#' L0
#' calculate the L0 norm of a vector
#' @param v the vector of which the L0 norm is calculated
#' @return the L0 norm of the vector v
#' @export L0


L0 <- function(v){
    sum(v != 0)
}

#' L2
#' calculate the L2 norm of a vector
#' @param v the vector of which the L2 norm is calculated
#' @return the L2 norm of the vector v
#' @export L2

L2 <- function(v){
    as.numeric(sqrt(sum(v*v)))
}


#' softThresh
#' returns the soft threshold of a vector
#' @param x the vector on which the soft threshold operation is performed
#' @param a the threshold value
#' @return the soft thresholding of vector x
#' @export

softThresh <- function(x,a){
    result <- numeric(length(x))
    result[which(x > a)] <- x[which(x > a)] - a
    result[which(x < -a)] <- x[which(x < -a)] + a
    return(result)
}



#' mst
#'
#' Build the adjancey matrix of an MST
#'
#'
#' @param X a distance matrix based on which the MST is built
#' @return the adjancey matrix of the MST
#' @examples
#' X1 <- runif(50, -1, 1)
#' X2 <- runif(50, -1, 1)
#' X <- cbind(X1, X2)
#' D <- dist(X)
#' mst(D)
#' @noRd


mst <- function (X){

    if (class(X) == "dist")
        X <- as.matrix(X)
    n <- dim(X)[1]
    N <- matrix(0, n, n)
    tree <- NULL
    large.value <- max(X) + 1
    diag(X) <- large.value
    index.i <- 1
    for (i in 1:(n - 1)) {
        tree <- c(tree, index.i)
        m <- apply(as.matrix(X[, tree]), 2, min)
        a <- sortIndex(X[, tree])[1, ]
        b <- sortIndex(m)[1]
        index.j <- tree[b]
        index.i <- a[b]
        N[index.i, index.j] <- 1
        N[index.j, index.i] <- 1
        for (j in tree) {
            X[index.i, j] <- large.value
            X[j, index.i] <- large.value
        }
    }
    dimnames(N) <- dimnames(X)
    class(N) <- "mst"
    return(N)


}


#' sortIndex
#'
#' sort each column of the given matrix X (default in the ascending order)
#'
#'
#' @param X a matrix/vector to be sorted
#' @param decrease sorted in ascending or descending order, default is FALSE
#' @return the sorted indices for each column vector
#' @examples
#'
#'
#' sortIndex(c(2,-1,0), decrease = TRUE)
#' sortIndex(cbind(c(2,3,1), c(5,10,-1)))
#'
#'
#' @export

sortIndex <- function(X, decrease = F){

    if (length(X) == 1){
        return(1)
    }

    if (!is.matrix(X)){
        X <- as.matrix(X)
    }

    indices <- apply(X, 2, function(v) order(rank(v)))

    if (decrease == T){
        indices <- apply(X, 2, function(v) order(rank(v), decreasing = T))
    }

    return(indices)

}


#' summary.hp
#'
#' OLS summary of the clusterwise regression (cwr) models
#'
#' @param object a cwr object returned by function HP
#' @param ... further arugments to be passed to
#' @return a list of "summary.lm" objects, with each object a summary of the linear model of a subgroup.
#'
#' @examples
#' n <- 50
#' p <- 3
#'
#' X <- matrix(rnorm(n * p), nrow = n)
#' Xj <- X[,1] # the threshold variable
#' beta1 <- rep(3,p)
#' beta2 <- rep(-3,p)
#'
#' index.g1 <- which(Xj <= 0)
#' index.g2 <- which(Xj > 0)
#'
#' y.g1 <- X[index.g1,] %*% beta1
#' y.g2 <- X[index.g2,] %*% beta2
#'
#' y <- rep(0,n)
#' y[index.g1] <- y.g1
#' y[index.g2] <- y.g2
#'
#' y <- y + rnorm(n = n, sd = 0.5)
#'
#' res <- HP(X, y, method = "mst", Xj = X[,1], max.no.cluster = 2)
#' s <- summary(res)
#'
#'  ## summary of the first subgroup model ##
#' s1 <- s$lm1.summary
#' coef(s1) # regression coefficient estimates of the first group
#' s1
#' @export

summary.hp <- function(object, ...){

    res.list <- object

    membership <- res.list$membership

    K <- length(res.list) - 1
    res.n <- table(membership)

    res <- list()

    for (k in 1:K){
        names(res.n)[k] <- paste("G", k, sep = "")
    }

    res$group.size <- res.n

    for (k in 1:K){
        res[[k+1]] <- summary.lm(res.list[[k+1]], ...)
        names(res)[k+1] <- paste("lm", k, ".summary", sep = "")
    }

    return(res)
}


summary.hp <- function(object, ...){

    res.list <- object

    membership <- res.list$membership

    K <- length(res.list) - 1
    res.n <- table(membership)

    res <- list()

    for (k in 1:K){
        names(res.n)[k] <- paste("G", k, sep = "")
    }

    res$group.size <- res.n

    for (k in 1:K){
        res[[k+1]] <- summary.lm(res.list[[k+1]], ...)
        names(res)[k+1] <- paste("lm", k, ".summary", sep = "")
    }

    return(res)
}






