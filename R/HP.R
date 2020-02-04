#' R.f
#'
#' primal proximal update of the PD algorithm
#'
#' @param x a vector of length p (for each obs)
#' @param y a scalar
#' @param beta a vector of length p (for each obs)
#' @param alpha step size
#' @param n sample size

#' @return updated result of the primal proximal
#' @noRd

R.f <- function(x, y, beta, alpha, n){

    p <- length(x)

    inverse <- solve(1/n * tcrossprod(x) +  1/alpha * diag(p))

    u <-  crossprod(t(inverse), (y/n * x  + 1/alpha * beta))

    return(u)
}


#' R.g
#'
#' dual proximal update of the PD algorithm
#'
#' @param u input vector
#' @return updated result of the dual proximal
#' @noRd

R.g.l2 <- function(u){
    min(1, 1/L2(u)) * u
}




#' iteration.cp
#'
#' Returns the individual beta estimates after penalization
#'
#' @param X the design matrix (the column vector of 1 should be included)
#' @param y the response vector
#' @param tau the tuning parameter
#' @param M the contrast matrix
#' @param B.init the initial value for individual beta's
#' @param maxit the maximum number of iterations; default = 1E5
#' @param epsilon the error tolerance; default = 1E-4
#' @return Returns the individual beta estimates after regularization
#' @examples
#'
#'
#' set.seed(121)
#'
#' #############
#' ## dataset ##
#' #############
#'
#' n <- 50
#' p <- 3
#' beta <- c(1,-1,1)
#' X <- matrix(rnorm(n * p), nrow = n)
#' e <- rnorm(n = n, sd = 0.5)
#' y <- X %*% beta + e
#' M <- build.M.c(n)
#'
#' iteration.cp(X, y, tau = sqrt(log(n)/n), M = M)
#' @importFrom Matrix sparseMatrix
#' @noRd

iteration.cp <- function(X, y,
                         tau = NULL, M,
                         B.init = NULL,
                         epsilon = 10^(-4),
                         maxit = 10^5){


    p <- ncol(X)
    n <- nrow(X)

    if (is.null(tau)){
        tau <- 0.1 * sqrt(log(n) / n)
    }

    if (is.null(B.init)){
        beta.ols <- coef(lm(y~X-1))
        B.init <- tcrossprod(rep(1,n), (beta.ols))
    }

    B.init <- t(B.init)


    # the step size
    M.sparse <- methods::as(M, "sparseMatrix")
    M.norm <- irlba::irlba(tau * M.sparse, nv = 1, maxit = 10^5)$d
    alpha <- 0.9 / M.norm

    ########################
    ## primal-dual update ##
    ########################


    M <- tau * M
    R.g <- R.g.l2

    U.init <- sign(tcrossprod(B.init, M))

    B.old <- B.init
    U.old <- U.init

    B.new <- matrix(0, nrow = p, ncol = n)

    term1 <- 1/ (2 * n) * L2( y - rowSums(X * t(B.new)) )
    term2 <- tau * sum(apply(tcrossprod(B.new, M), 2, L2) )
    term.old <-term1 + term2


    for (j in 1:maxit){

        # step 1:

        # matrix of dim = p * n
        B.input <- B.old - alpha * tcrossprod(U.old, t(M))

        # matrix of dim = p * n
        for (i in 1:n){
            B.new[,i] <- R.f(x = X[i,], y = y[i],
                             beta = B.input[,i], n = n,
                             alpha = alpha)
        }

        # step 2:
        # matrix of dim = p * S
        U.input <- U.old + alpha * tcrossprod( (2 * B.new - B.old), M)


        #matrix of dim = p * S
        U.new <- apply(U.input, 2, R.g)


        term1 <- 1/ (2 * n) * L2( y - rowSums(X * t(B.new)) )
        term2 <- tau * sum(apply(tcrossprod(B.new, M), 2, L2) )

        term.new <- term1 + term2

        if (abs((term.new - term.old)/ term.old) < epsilon){
            break
        }

        # update:
        B.old <- B.new
        U.old <- U.new
        term.old <- term.new
    }

    return(t(B.new))
}



#' cluster
#'
#' step 1 update of C-P algorithm
#'
#' @param B.hat the matrix of case-specific beta's
#' @param C number of clusters
#' @return membership labels
#' @noRd


cluster <- function(B.hat, C){

    d <- stats::dist(B.hat, method = "euclidean")
    hc <- stats::hclust(d, method = "complete")

    membership <- stats::cutree(hc, k = C)
    return(membership)
}



#' HP
#'
#' Perform heterogeneity pursuit on the given dataset
#'
#' @param X the design matrix; intercept column should bes included, if any.
#' @param y the response vector of length n.
#' @param tau the tuning parameter of the penalized regression; if not provided, tau = 0.1 * sqrt(n / log(n)).
#' @param method method = c(``latent'', ``mst'', ``knn'').
#' @param Xj the threshold variable, must be provided if method = ``mst'' or ``knn''
#' @param K the number of neighbors in a KNN graph; It may be specfied if method = "knn" or "latent"; default value is 3.
#' @param  no.cluster.search TRUE or FALSE; if TRUE, the optimal number of clusters is decided based on some information criterion.
#' @param IC  IC = c("aic", "bic", "mdl"), the information criterion based on which the cluster number of the dataset is determined; default is IC = "bic".
#' @param max.no.cluster If no.cluster.search = TRUE, it is the max number of clusters which the data is partitioned into; if no.cluster.search = FALSE, max.no.cluster is the user-specified cluster numbers. The default value is max(n^(1/3), 5).
#' @return a list that contains (i). a vector of membership indicators; (ii). the OLS fit of each sub-model.
#' @examples
#'
#'
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
#' HP(X, y, method = "mst", Xj = X[,1], max.no.cluster = 2)
#' HP(X, y, method = "knn", Xj = X[,1], max.no.cluster = 2)
#' HP(X, y, method = "latent", max.no.cluster = 2)
#' HP(X, y, method = "knn", Xj = X[,1])
#' @export HP

HP <- function(X, y, tau = NULL, method,
               Xj = NULL, K = 3, no.cluster.search = FALSE,
               IC = "bic",
               max.no.cluster = NULL){

    if (method == "knn" & (is.null(Xj)) ){
        stop("please specify the threshold variable.")
    }

    if (method == "mst" & (is.null(Xj)) ){
        stop("please specify the threshold variable.")
    }

    n <- nrow(X)
    p <- ncol(X) + 1

    cl <- max.no.cluster

    if (is.null(cl)){

        cl <- floor(min(n^(1/3), 5))
    }

    if (method == "mst" & (!is.null(Xj))){
        M <- build.M.mst(Xj)
    }

    if (method == "knn" & (!is.null(Xj))){
        D <- as.matrix(dist(Xj))
        M <- build.M.nbr(K, D = D)
    }

    if (no.cluster.search){

        membership <- matrix(0, ncol = n, nrow = cl)
        ic <- rep(0, cl)

        for (l in 1:cl){

            if (method == "latent"){
                nbr.K <- sch(X, y, rs = 30, k  = l)
                M <- build.M.nbr(K, nbr.K)
            }

            B.hat <- iteration.cp(X = X, y = y, tau = tau, M = M)

            membership[l,] <- cluster(B.hat, l)
            ic[l] <- ic(X, y, membership[l,], type = IC)
        }

        index.chosen <- which.min(ic)
        membership.chosen <- membership[index.chosen,]

    }

    if (!no.cluster.search){

        if (method == "latent"){

            nbr.K <- sch(X, y, rs = 30, k  = cl)
            M <- build.M.nbr(K, nbr.K)

        }

        B.hat <- iteration.cp(X = X, y = y, tau = tau, M = M)

        membership.chosen <- cluster(B.hat, cl)
    }

    cc <- length(unique(membership.chosen))

    res.list <- list(membership = membership.chosen)
    model <- list()

    for (c in 1:cc){

        index.G <- which(membership.chosen == c)

        if(length(index.G) < p){
            warning(sprintf("sub-group size too small to fit OLS."))
            return(res.list)
        }

        y.G <- y[index.G]
        X.G <- X[index.G,]
        model[[(c)]] <- lm(y.G ~ X.G - 1)
        names(model)[c] <- paste("lm", c, sep = "")
    }

    res.list$model <- model
    class(res.list) <- "cwr"

    return(res.list)
}






