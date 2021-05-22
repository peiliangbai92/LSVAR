#' Function to generate a VAR process
#' @description A function to generate synthetic time series process based on the given structure
#' @param n the length of time series
#' @param p the number of multivariate time series
#' @param struct a character string indicating the structure of the transition matrix, here are three options:
#' sparse, low rank and LS (low rank plus sparse)
#' @param sp_density a numeric value, indicating the sparsity density of sparse components, default is 0.1
#' @param signal a numeric value, indicating the magnitude of transition matrix
#' @param rank a positive integer, the rank for low rank component
#' @param singular_vals a numeric vector, indicating the singular values for the low rank component, the length of
#' singular value must equal to the rank
#' @param spectral_radius a numeric value, controlling the stability of the process, default is 0.9
#' @param sigma a numeric matrix, indicating the covariance matrix of noise term
#' @param skip a numeric value, indicating the number of skipped time points in the beginning of the process
#' @param seed an integer, indicating the seed for random seed.
#' @import mvtnorm
#' @import pracma
#' @importFrom igraph erdos.renyi.game
#' @importFrom igraph get.adjacency
#' @importFrom igraph V
#' @importFrom stats rnorm
#' @return A list object, including
#' \describe{
#'     \item{series}{the generated time series}
#'     \item{noise}{the noise term}
#'     \item{model_param}{true transition matrix}
#' }
#' @export
#' @examples
#' n <- 300; p <- 15
#' signal <- 0.75
#' rank <- 3
#' singular_vals <- c(1, 0.75, 0.5)
#' try <- testVAR(n, p, struct = "LS", signal = signal, rank = rank,
#'                singular_vals = singular_vals)
#' data <- as.matrix(try$series)
testVAR <- function(n, p, struct = c("sparse", "low rank", "LS")[1], sp_density = 0.1,
                    signal = NULL, rank = NULL, singular_vals, spectral_radius = 0.9,
                    sigma = NULL, skip = 50, seed = 1){
    if(!(struct %in% c("sparse", "lowrank", "LS"))){
        stop("Error in strcture, must be in sparse, lowrank, LS!")
    }

    if((struct == "lowrank" | struct == "LS") & is.null(rank)){
        stop("Error! Should provide rank for the low rank component!")
    }

    if(spectral_radius >= 1){
        stop("Error! The spectral radius should be less than 1!")
    }

    if(is.null(sigma)){
        sigma <- 0.01*diag(p)
    }

    if(!is.matrix(sigma)){
        sigma <- as.matrix(sigma)
    }


    if(struct == "sparse"){
        if(is.null(signal)){
            stop("Error: should provide signal for sparse component!")
        }
        set.seed(seed)
        g <- erdos.renyi.game(p, round(sp_density * p^2), type = "gnm", directed = TRUE)
        adj_mat <- as.matrix(get.adjacency(g))
        sparse_mat <- adj_mat * signal

        # check stationarity
        max_eigen <- max(abs(eigen(sparse_mat)$values))
        cat(paste("Spectral radius is:", max_eigen, "\n", sep = " "))
        if(max_eigen >= 1){
            sparse_mat <- sparse_mat * spectral_radius / max_eigen
        }
        phi <- sparse_mat
    }

    if(struct == "low rank"){
        set.seed(seed)
        L_basis <- randortho(p)
        lowrank_mat <- matrix(0, p, p)
        for(rk in 1:rank){
            lowrank_mat <- lowrank_mat + singular_vals[rk] * (L_basis[,rk] %*% t(L_basis[,rk]))
        }

        # check stationarity
        max_eigen <- max(abs(eigen(lowrank_mat)$values))
        cat(paste("Spectral radius is:", max_eigen, "\n", sep = " "))
        if(max_eigen >= 1){
            lowrank_mat <- lowrank_mat * spectral_radius / max_eigen
        }
        phi <- lowrank_mat
    }

    if(struct == "LS"){
        set.seed(seed)
        L_basis <- randortho(p)
        lr_comp <- matrix(0, p, p)
        for(rk in 1:rank){
            lr_comp <- lr_comp + singular_vals[rk] * (L_basis[,rk] %*% t(L_basis[,rk]))
        }
        g <- erdos.renyi.game(p, round(sp_density * p^2), type = "gnm", directed = TRUE)
        adj_mat <- as.matrix(get.adjacency(g))
        sp_comp <- adj_mat * signal
        phi <- lr_comp + sp_comp

        # check stationarity
        max_eigen <- max(abs(eigen(phi)$values))
        cat(paste("Spectral radius is:", max_eigen, "\n", sep = " "))
        if(max_eigen >= 1){
            phi <- phi * spectral_radius / max_eigen
        }
    }

    # generate process
    N <- n + skip
    series <- matrix(0, nrow = N, ncol = p)
    noise <- rmvnorm(N, rep(0, p), sigma)
    series[1,] <- rnorm(p, 0, 1)
    for(t in 2:N){
        series[t, ] <- series[t-1, ] %*% phi + noise[t, ]
    }
    series <- series[-(1:skip),]
    noise <- noise[-(1:skip),]

    return(list(series = series, noise = noise, model_param = phi))
}



#' A function to solve low rank plus sparse model estimation using FISTA algorithm
#' @description A function to solve low rank plus sparse model estimation
#' @param data A numeric dataset with size of n by p
#' @param lambda A positive numeric value, indicating the tuning parameter for sparse component
#' @param mu A positive numeric value, indicating the tuning parameter for low rank component
#' @param alpha_L The constraint coefficient of low rank component, default is 0.25
#' @param niter The maximum number of iterations required for FISTA
#' @param backtracking A boolean argument, indicating that use backtracking in the FISTA
#' @param x.true A p by p matrix, the true model parameter. Only available for simulation.
#' @return A S3 object of class \code{LSVAR}, including
#' \describe{
#'   \item{est_phi}{estimated model parameter}
#'   \item{sparse.comp}{Estimated sparse component}
#'   \item{lr.comp}{Estimated low-rank component}
#'   \item{obj.val}{Values of objective function}
#'   \item{rel.err}{Relative errors compared with the true model parameters if available}
#' }
#' @export
#' @examples
#' n <- 300
#' p <- 20
#' try <- testVAR(n, p, struct = "LS", signal = 0.75, rank = 2,
#'                singular_vals = c(1, 0.8))
#' data <- as.matrix(try$series)
#' lambda <- 0.1; mu <- 1
#' fit <- fista.LpS(data, lambda = lambda, mu = mu, x.true = try$model_param)
#' summary(fit, threshold = 0.2)
fista.LpS <- function(data, lambda, mu, alpha_L = 0.25,
                      niter = 100, backtracking = TRUE, x.true){
    n <- dim(data)[1]
    p <- dim(data)[2]
    A <- data[1:(n-1), ]
    b <- data[2:n, ]
    tnew = t <- 1
    x1 <- matrix(0, nrow = p, ncol = p)
    xnew1 = xnew2 <- x1
    y1 = y2 <- x1
    AtA <- t(A) %*% A
    Atb <- t(A) %*% b

    if(backtracking == TRUE){
        L <- norm(A, "F")^2 / 5
        gamma <- 2
    }else{
        L <- norm(A, "F")^2
    }

    obj.val = rel.err <- c()
    for(i in 1:niter){
        if(backtracking == TRUE){
            L.bar <- L
            found <- FALSE
            while(found == FALSE){
                y <- y1 + y2
                prox1 <- prox.sparse.func(y1, y, A, b, 2*L.bar, lambda, AtA, Atb)
                prox2 <- prox.nuclear.func(y2, y, A, b, 2*L.bar, mu, AtA, Atb)

                ### Restricted solution space
                for(j in 1:p){
                    for(k in 1:p){
                        if(abs(prox2[j,k]) > alpha_L){
                            prox2[j,k] <- alpha_L * sign(prox2[j,k])
                        }
                    }
                }

                prox <- prox1 + prox2
                if(f.func(prox, A, b) <= Q.func(prox, y, A, b, L.bar, AtA, Atb)){
                    found <- TRUE
                }else{
                    L.bar <- L.bar * gamma
                }
            }
            L <- L.bar
        }
        x1 <- xnew1
        x2 <- xnew2
        xnew1 <- prox1
        xnew2 <- prox2
        t = tnew
        tnew <- (1 + sqrt(1 + 4*t^2))/2
        y1 <- xnew1 + (t - 1) / tnew * (xnew1 - x1)
        y2 <- xnew2 + (t - 1) / tnew * (xnew2 - x2)
        xnew <- xnew1 + xnew2

        obj.val <- c(obj.val, f.func(xnew, A, b) + sparse.pen(xnew1, lambda) + nuclear.pen(xnew2, mu))
        rel.err <- c(rel.err, norm(xnew - x.true, "F") / norm(x.true, "F"))
    }
    final.result <- structure(list(est_phi = xnew1+xnew2,
                                   sparse.comp = xnew1,
                                   lr.comp = xnew2,
                                   obj.val = obj.val,
                                   rel.err = rel.err), class = "LSVAR")
    return(final.result)
}


#' plot sparse component for use igraph and network layout
#' @description Plot a network to illustrate the estimated sparse component
#' @param mat a p by p matrix, indicating the sparse component
#' @param threshold the threshold for presenting the edges in the network
#' @importFrom igraph plot.igraph
#' @importFrom igraph graph.adjacency
#' @importFrom igraph layout_in_circle
#' @importFrom igraph plot.igraph
#' @importFrom igraph V
#' @return A network plot for the sparse component
#' @export
#' @examples
#' set.seed(1)
#' est_mats <- matrix(rnorm(400, 0, 1), 20, 20)
#' plot_network(est_mats, threshold = 0.1)
plot_network <- function(mat, threshold = 0.1){
    p <- dim(mat)[1]
    adj_mat <- matrix(0, p, p)
    for(r in 1:p){
        for(c in 1:p){
            if(abs(mat[r,c]) > threshold){
                adj_mat[r,c] <- 1
            }
        }
    }
    net <- graph.adjacency(adj_mat, "directed", diag = FALSE)
    l <- layout_in_circle(net)
    plot.igraph(net, vertex.label = V(net)$name, layout = l, vertex.label.font = 2,
                vertex.label.color = "black", edge.color = "gray40", edge.arrow.size = 0.1,
                vertex.shape ="circle", vertex.color = "light blue", vertex.label.cex = 0.8)
}


#' Summary of LSVAR S3 class
#' @description summary function for S3 class for the fitting result
#' @param object the S3 class object of \code{LSVAR}
#' @param threshold the threshold for sparse component visualization
#' @param ... not in use
#' @return A series of summary for the S3 result
#' @export
#' @examples
#' n <- 300
#' p <- 20
#' try <- testVAR(n, p, struct = "LS", signal = 0.75, rank = 2,
#'                singular_vals = c(1, 0.8))
#' data <- as.matrix(try$series)
#' lambda <- 0.1; mu <- 1
#' fit <- fista.LpS(data, lambda = lambda, mu = mu, x.true = try$model_param)
#' summary(fit, threshold = 0.2)
summary.LSVAR <- function(object, threshold = 0.2, ...){
    est_phi <- object$est_phi
    sparse <- object$sparse.comp
    lowrank <- object$lr.comp
    re <- object$rel.err
    p <- dim(est_phi)[1]
    sp_mat <- matrix(0, p, p)
    for(r in 1:p){
        for(c in 1:p){
            if(abs(sparse[r,c]) > threshold){
                sp_mat[r,c] <- sparse[r,c]
            }
        }
    }
    cat("===========================================\n")
    cat(paste("Rank for the estimated low rank component:", qr(lowrank)$rank, "\n", sep = " "))
    cat("Print the relative error curve: \n")
    cat("Plot the estimated sparse component: \n")
    plot_network(sparse, threshold = threshold)
    plot(re, type = 'l')
    cat("===========================================\n")
}




#' Shrinkage function for sparse soft-thresholding
#' @description Shrinkage function for sparse soft-thresholding
#' @param y A matrix, or a vector for thresholding
#' @param tau A positive number, threshold
#' @return A thresholded matrix, or vector
shrinkage <- function(y, tau){
    z <- matrix(0, nrow = nrow(y), ncol = ncol(y))
    for(i in 1:nrow(y)){
        for(j in 1:ncol(y)){
            z[i,j] <- sign(y[i,j]) * max(abs(y[i,j]) - tau, 0)
        }
    }
    return(z)
}

#' Shrinkage function for low-rank soft-thresholding
#' @description Shrinkage function for low-rank soft-thresholding
#' @param y A matrix, or a vector for thresholding
#' @param tau A positive number, threshold
#' @return A thresholded matrix, or vector
shrinkage.lr <- function(y, tau){
    z <- rep(0, length(y))
    for(i in 1:length(y)){
        z[i] <- sign(y[i]) * max(0, abs(y[i]) - tau)
    }
    return(z)
}


#' Gradient function of quardratic loss
#' @description Gradient function of quardratic loss
#' @param x A vector, or matrix, indicating the model parameter
#' @param AtA A p by p Gram matrix for corresponding design matrix A
#' @param Atb An inner product for design matrix A and corresponding matrix (vector) b
#' @return Value of gradients
gradf.func <- function(x, AtA, Atb){
    return(AtA %*% x - Atb)
}


#' Nuclear norm penalty for low-rank component
#' @description Nuclear norm penalty for low-rank component
#' @param x Model parameter
#' @param lambda Tuning parameter
#' @return Value of nuclear norm penalty term
nuclear.pen <- function(x, lambda){
    d <- svd(x)$d
    return(lambda * sum(d))
}


#' L1-norm penalty for sparse component
#' @description L1-norm penalty for sparse component
#' @param x Model parameter
#' @param lambda Tuning parameter
#' @return Value of l1-norm penalty term
sparse.pen <- function(x, lambda){
    return(lambda*sum(x))
}


#' Main loss function for quardratic loss
#' @description Main loss function
#' @param x Model parameters
#' @param A Design matrix with size of n by p
#' @param b Correspond vector or matrix
#' @return Value of objective function
f.func <- function(x, A, b){
    return(0.5 * norm(A %*% x - b, "F")^2)
}


#' An auxiliary function in FISTA algorithm
#' @description Auxiliary function for FISTA implementation
#' @param x Model parameter for previous update
#' @param y Model parameter for updating
#' @param A An n by p design matrix
#' @param b A correspond vector, or matrix with size of n by 1 or n by p
#' @param L Learning rate
#' @param AtA Gram matrix for design matrix A
#' @param Atb Inner product for design matrix A and correspond vector b
#' @return Value of function Q
Q.func <- function(x, y, A, b, L, AtA, Atb){
    return(f.func(y, A, b) + sum((x - y) * gradf.func(y, AtA, Atb)) + 0.5 * L * norm(x - y, "F")^2)
}


#' Proximal function with nuclear norm penalty updating
#' @description Proximal function with nuclear norm
#' @param w1 previously updated model parameter
#' @param y updated model parameter
#' @param A design matrix
#' @param b correspond vector, or matrix
#' @param L learning rate
#' @param lambda tuning parameter for low-rank component
#' @param AtA Gram matrix of design matrix A
#' @param Atb inner product of design matrix A and correspond vector b
#' @return Value of proximal function with nuclear norm penalty
prox.nuclear.func <- function(w1, y, A, b, L, lambda, AtA, Atb){
    Y <- w1 - (1 / L) * gradf.func(y, AtA, Atb)
    d <- shrinkage.lr(svd(Y)$d, 2*lambda / L)
    return(svd(Y)$u %*% diag(d) %*% t(svd(Y)$v))
}


#' Proximal function with l1-norm penalty updating
#' @description Proximal function with l1-norm
#' @param w1 previously updated model parameter
#' @param y updated model parameter
#' @param A design matrix
#' @param b correspond vector, or matrix
#' @param L learning rate
#' @param lambda tuning parameter for sparse component
#' @param AtA Gram matrix of design matrix A
#' @param Atb inner product of design matrix A and correspond vector b
#' @return Value of proximal function with l1-norm penalty
prox.sparse.func <- function(w1, y, A, b, L, lambda, AtA, Atb){
    Y <- w1 - (1 / L) * gradf.func(y, AtA, Atb)
    return(shrinkage(Y, 2*lambda / L))
}



#' Objective function
#' @description objective function, main loss function and penalties
#' @param x.lr low-rank component
#' @param x.sparse sparse component
#' @param A design matrix
#' @param b correspond vector
#' @param lambda a tuning parameter for sparse component
#' @param mu a tuning parameter for low-rank component
#' @return value of objective function
obj.func <- function(x.lr, x.sparse, A, b, lambda, mu){
    ### x.sparse is a list
    m <- length(x.sparse)
    loss <- 0
    for(i in 1:m){
        loss <- loss + f.func((x.lr[[i]] + x.sparse[[i]]), A, b) + sparse.pen(x.sparse[[i]], lambda) + nuclear.pen(x.lr[[i]], mu)
    }
    return(loss)
}
