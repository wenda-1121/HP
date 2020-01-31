#' build.M.c
#'
#' Build the contrast matrix based on a complete graph with n nodes
#'
#'
#' @param n number of nodes for the complete graph
#' @return Returns the contrast matrix of a complete graph with n nodes
#' @examples
#'
#'
#' build.M.c(10)
#'
#' @import stats
#' @export build.M.c


build.M.c <- function(n){

    M.mini <-  matrix(rep(0, n*choose(n,2)), nrow = choose(n,2))

    k <- 1
    for (i in 1:(n-1)){
        for (j in (i+1):n){
            m.mini <- rep(0, n)
            m.mini[i] <- 1
            m.mini[j] <- -1
            M.mini[k,] <- m.mini
            k <- k+1
        }
    }

    return(M.mini)
}



#' build.M.nbr
#'
#' Build the contrast matrix of an KNN graph with n nodes
#'
#'
#' @param K number of neighbors for each node
#' @param nbr.K nbr.K[i,] = the node indices ranked by proximity to node i
#' (from the nearest to the farthest)
#' @param D the distance matrix
#' @return Returns the contrast matrix of a KNN graph with n nodes
#' @examples
#'
#'
#' set.seed(121)
#'
#' n <- 5 # number of nodes
#' X <- runif(n, -1, 1) # nodes are generated from U(-1, 1)
#' D <- as.matrix(dist(X)) # distance matrix for the nodes
#' nbr.K <- t(sortIndex(D)) # ranked indices based on distance
#' M1 <- build.M.nbr(K = 2, nbr.K = nbr.K[,-1])
#' M2 <- build.M.nbr(K = 2, D = D)
#'
#' mean(M1 == M2) == 1 #returns TRUE
#'
#' @export build.M.nbr


build.M.nbr <- function(K, nbr.K = NULL, D = NULL){

    if(!is.null(D)){

        D <- as.matrix(D)
        nbr.K <- t(sortIndex(D))
        nbr.K <- nbr.K[,-1]
    }

    n <- nrow(nbr.K)

    M.mini <- matrix(rep(0,n*K*n), ncol = n)

    for (i in 1:n){

        start <- (i-1) * K + 1
        end <- i * K
        r <- 1

        for (l in start:end){
            M.mini[l,i] <- 1
            j <- nbr.K[i,r]
            M.mini[l,j] <- -1
            r <- r + 1
        }
    }

    # repetitive edges
    index.rep <- NULL

    for (i in 1:nrow(M.mini)){

        Mi <- M.mini[i,]

        for (j in i:nrow(M.mini)){
            Mj <- M.mini[j,]

            if (mean(Mj == -Mi) == 1){
                index.rep <- c(index.rep, j)
            }
        }
    }

    M.mini <- M.mini[-index.rep,]

   return(M.mini)
}


#' build.M.mst
#'
#' Build the contrast matrix of an MST constructed based on the threshold
#' variable provided.
#'
#'
#' @param Xj the threshold variable
#' @return Returns the contrast matrix of a MST constructed based on the
#' threshold variable provided
#' @examples
#'
#' set.seed(121)
#' n <- 20
#' Xj <- runif(n, -1 ,1)
#'
#' build.M.mst(Xj)
#'
#' @export build.M.mst
#'


build.M.mst <- function(Xj){

    Xj <- as.matrix(Xj)
    n <- nrow(Xj)
    d <- as.matrix(dist(Xj))
    MST <- mst(d)

    MST[lower.tri(MST, diag = FALSE)] <- 0

    M.mst <- matrix(0, nrow = n-1, ncol = n)

    k <- 1

    for (i in 1:n){

        # if not all entries in the ith row are 0
        if ( mean(MST[i,] == 0) != 1){

            node.index <- which(MST[i,] !=0)
            l <- length(node.index)
            i.index <- i

            for (j in 1:l){
                j.index <- node.index[j]
                M.mst[k,i.index] <- 1
                M.mst[k,j.index] <- -1
                k <- k + 1
            }

        }
    }

    return(M.mst)
}





