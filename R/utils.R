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


#' summary.cwr
#'
#' OLS summary of the clusterwise regression (cwr) models
#'
#' @param model.list a list of lm objects
#' @return a list of length K, where K is the number of clusters and each list element is an OLS summary of a subgroup model
#'
#' @export

summary.cwr <- function(model.list){

    K <- length(model.list)

    res <- list()

    for (k in 1:K){
        res[[k]] <- summary(model.list[[k]])
        names(res)[k] <- paste("lm", k, ".summary", sep = "")
    }

    return(res)
}









