#https://www.r-bloggers.com/2019/07/sampling-paths-from-a-gaussian-process/

library(MASS)
# generate covariance matrix for points in `x` using given kernel function
cov_matrix <- function(x, kernel_fn, ...) {
  outer(x, x, function(a, b) kernel_fn(a, b, ...))
}


# given x coordinates, take N draws from kernel function at those points
draw_samples <- function(x, N, seed = 123, kernel_fn, centred=FALSE,squared=FALSE, ...) {
  Y <- matrix(NA, nrow = length(x), ncol = N)
  C<-diag(length(x))-1/length(x)*matrix(1,nrow=length(x),ncol = length(x))
  set.seed(seed)
  for (n in 1:N) {
    K <- cov_matrix(x, kernel_fn, ...)
    if (centred==TRUE & squared==FALSE){
      K <- C%*%K%*%C
    } else if (centred==FALSE & squared==TRUE){
      K <- K%*%K
    } else if (centred==TRUE & squared==TRUE){
      K <- C%*%K%*%C
      K <- K%*%K
    }
    Y[, n] <- mvrnorm(1, mu = rep(0, length(x)), Sigma = K)
  }
  Y
}


cov_matrix_3D <- function(x,y, kernel_fn, ...) {
  outer(x, y, function(a, b) kernel_fn(a, b, ...))
}


draw_samples_3D <- function(x1a,x2a, N, seed = 123, kernel_fn, centred=FALSE,squared=FALSE, ...) {
  X <- expand.grid(x1a, x2a)
  x1<-X$Var1
  x2<-X$Var2
  Y <- matrix(NA, nrow=N, ncol = dim(X)[1])
  C<-diag(length(x1))-1/length(x1)*matrix(1,nrow=length(x1),ncol = length(x1))
  set.seed(seed)
  K1 <- cov_matrix_3D(x1,x1, kernel_fn, ...)
  K2 <- cov_matrix_3D(x2,x2, kernel_fn, ...)
  if (centred==TRUE & squared==FALSE){
    K1 <- C%*%K1%*%C
    K2 <- C%*%K2%*%C
  } else if (centred==FALSE & squared==TRUE){
    K1 <- K1%*%K1
    K2 <- K2%*%K2
  } else if (centred==TRUE & squared==TRUE){
    K1 <- C%*%K1%*%C
    K2 <- C%*%K2%*%C
    K1 <- K1%*%K1
    K2 <- K2%*%K2
  }

  Y <- mvrnorm(N, mu=rep(0,length(x1)),Sigma = (1*(matrix(1,nrow=length(x1),ncol=length(x1))+K1+K2)))
  return(Y)
}

# Kernels
# Browian Motion

bm_kernel <- function(x, y) {
  1/2*(abs(x)+abs(y)-abs(x-y))
}


# rbf kernel
se_kernel <- function(x, y, sigma = 1, length = 1) {
  sigma^2 * exp(- (x - y)^2 / (2 * length^2))
}

# Matern
matern_kernel <- function(x, y, nu = 1.5, sigma = 1, l = 1) {
  if (!(nu %in% c(0.5, 1.5, 2.5))) {
    stop("p must be equal to 0.5, 1.5 or 2.5")
  }
  p <- nu - 0.5
  d <- abs(x - y)
  if (p == 0) {
    sigma^2 * exp(- d / l)
  } else if (p == 1) {
    sigma^2 * (1 + sqrt(3)*d/l) * exp(- sqrt(3)*d/l)
  } else {
    sigma^2 * (1 + sqrt(5)*d/l + 5*d^2 / (3*l^2)) * exp(-sqrt(5)*d/l)
  }
}


# k <- function(x,y=NULL, kernel_fn, centred=FALSE,squared=FALSE,...) {
#   if (is.null(y)){
#     y<-x
#   }
#   n <- nrow(x)
#   m <- nrow(y)
#   K <- matrix(0, n, m)
#     for (i in 1:n) {
#       for (j in 1:m) {
#         K[i, j] <- kernel_fn(x[i, ], x[j, ], ...)
#       }
#     }
#   C<-diag(n)-1/n*matrix(1,nrow=n,ncol = n)
#   if (centred==TRUE & squared==FALSE){
#     K <- t(C)%*%K%*%C
#   } else if (centred==FALSE & squared==TRUE){
#     K <- t(K)%*%K
#   } else if (centred==TRUE & squared==TRUE){
#     K <- t(C)%*%K%*%C
#     K <- t(K)%*%K
#   }
#   K
# }
# 
# 
# 
# 
# gram_matrix <- function(x, kernel_fn, centred=FALSE,squared=FALSE,...) {
#   n <- nrow(x)
#   if (is.vector(x)){
#     K <- outer(x, x, function(a, b) kernel_fn(a, b, ...))
#   } else if (is.matrix(x)) {
#     K <- matrix(0, n, n)
#     for (i in 1:n) {
#       for (j in 1:n) {
#         K[i, j] <- kernel_fn(x[i, ], x[j, ], ...)
#       }
#     }
#   }
#   C<-diag(n)-1/n*matrix(1,nrow=n,ncol = n)
#   if (centred==TRUE & squared==FALSE){
#     K <- t(C)%*%K%*%C
#   } else if (centred==FALSE & squared==TRUE){
#     K <- t(K)%*%K
#   } else if (centred==TRUE & squared==TRUE){
#     K <- t(C)%*%K%*%C
#     K <- t(K)%*%K
#   }
#   K
# }




#See https://github.com/haziqj/iprior/blob/master/R/Plots.R line 268
# p <- ggplot() +
#   scale_x_continuous(breaks = NULL, name = expression(italic(y))) +
#   scale_y_continuous(breaks = NULL) +
#   geom_line(data = melted.ppc, stat = "density", alpha = 0.5,
#             aes(x = value, group = variable, col = "yrep", size = "yrep")) +
#   geom_line(aes(x = get_y(x), col = "y", size = "y"), stat = "density") +
#   theme(legend.position = "bottom") +
#   scale_colour_manual(
#     name = NULL, labels = c("Observed", "Replications"),
#     values = c("grey10", "steelblue3")
#   ) +
#   scale_size_manual(
#     name = NULL, labels = c("Observed", "Replications"),
#     values = c(1.1, 0.19)
#   ) +
#   labs(y = "Density", title = "Posterior predictive density check") +
#   theme_bw() +
#   theme(legend.position = c(0.9, 0.5))

#col_list <- c("red", "blue", "black")  # for line colors


# centered squared kernel but takes long computationally
# bm_kernel_centred_squared <- function(x, y) {
#   n_x <- length(x)
#   n_y <- length(y)
#   k1 <- sum(sapply(x, function(i) bm_kernel(x, i)))
#   k2 <- sum(sapply(y, function(j) bm_kernel(y, j)))
#   k3 <- sum(outer(x, y, bm_kernel))
#   
#   result <- (bm_kernel(x, y) - (1 / n_x) * k1 - (1 / n_y) * k2 + (1 / (n_x * n_y)) * k3)^2
#   return(result)
# }

