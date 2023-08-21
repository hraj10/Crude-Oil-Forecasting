#https://www.r-bloggers.com/2019/07/sampling-paths-from-a-gaussian-process/
library(vscDebugger)
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
    if (centred==TRUE && squared==FALSE){
      K <- C%*%K%*%C
    } else if (centred==FALSE && squared==TRUE){
      K <- K%*%K
    } else if (centred==TRUE && squared==TRUE){
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
  if (centred==TRUE && squared==FALSE){
    K1 <- C%*%K1%*%C
    K2 <- C%*%K2%*%C
  } else if (centred==FALSE && squared==TRUE){
    K1 <- K1%*%K1
    K2 <- K2%*%K2
  } else if (centred==TRUE && squared==TRUE){
    K1 <- C%*%K1%*%C
    K2 <- C%*%K2%*%C
    K1 <- K1%*%K1
    K2 <- K2%*%K2
  }
  
  Y <- mvrnorm(N, mu=rep(0,length(x1)),Sigma = (matrix(1,nrow=length(x1),ncol=length(x1))+K1+K2))
  return(Y)
}

# Kernels
# Browian Motion

bm_kernel <- function(x, y) {
    1/2*(abs(x)+abs(y)-abs(x-y))
}
bm_kernel_L2 <- function(x,y){
  0.5 * (sqrt(sum(abs(x)^2)) + sqrt(sum(abs(y)^2)) - sqrt(sum(abs(x - y)^2)))
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




#a<-c(1,2,3,4)
#b<-cbind(a,b)
#c<-cbind(a,rep(0,4))
#gram_matrix(b,bm_kernel_L2, TRUE,FALSE) 
#library(iprior)
#kern_fbm(b)

# k <- function(x,y=NULL, kernel_fn, centred=FALSE,squared=FALSE,...) {
#   if (is.null(y)){
#     y<-x
#   }
#   x<-matrix(x)
#   y<matrix(y)
#   n <- nrow(x)
#   m <- nrow(y)
#   K <- matrix(0, n, m)
#   for (i in 1:n) {
#     for (j in 1:m) {
#       K[i, j] <- kernel_fn(x[i, ], x[j, ], ...)
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
#k(a, y=NULL, bm_kernel, TRUE, FALSE)
#kern_fbm(a)

K_gram <- function(x,y=NULL, kernel_fn, centred=FALSE,squared=FALSE,...) {
  x<-as.matrix(x)
  if (is.null(y)){
    y<-x
  } else {
    y<-as.matrix(y)
  }
  n <- nrow(x)
  m <- nrow(y)
  K <- matrix(0, n, m)
  for (i in 1:n) {
    for (j in 1:m) {
        K[i, j] <- kernel_fn(x[i, ], y[j, ], ...)
      }
  }
  C<-diag(n)-1/n*matrix(1,nrow=n,ncol = n)
  if (centred==TRUE && squared==FALSE){
    K <- t(C)%*%K%*%C
  } else if (centred==FALSE && squared==TRUE){
    K <- t(K)%*%K
  } else if (centred==TRUE && squared==TRUE){
    K <- t(C)%*%K%*%C
    K <- t(K)%*%K
  }
  K
}
#K_gram(as.matrix(a),y=as.matrix(a), bm_kernel_L2, centred=TRUE, squared=FALSE)

k_plain <- function(x,y=NULL, kernel_fn,...) {
  x<-as.matrix(x)
  if (is.null(y)){
    y<-x
  } else {
    y<-as.matrix(y)
  }
  n <- nrow(x)
  m <- nrow(y)
  K <- matrix(0, n, m)
  for (i in 1:n) {
    for (j in 1:m) {
      K[i, j] <- kernel_fn(x[i, ], y[j, ], ...)
    }
  }
  K
}

# kernel_vec<- function(x,y, kernel_fn, squared=FALSE,...){
#   n<-nrow(as.matrix(x))
#   k_orig<-k_plain(x,y,kernel_fn = bm_kernel_L2)
#   print(k_orig)
#   k<-k_orig
#   K_gram_matrix<-k_plain(x,y=NULL, bm_kernel_L2)
#   first_comp<-as.matrix(colSums(K_gram_matrix))
#   print(first_comp)
#   second_comp<-colSums(k_orig)
#   print(second_comp)
#   third_comp<-sum(K_gram_matrix)
#   print(third_comp)
#   cent<-k_orig-1/n*first_comp-1/n*second_comp+1/(n^2)*third_comp
#   if (squared==TRUE) {
#     cent=K_gram(x,y=NULL, bm_kernel_L2, centred=TRUE, squared=FALSE)%*%cent
#   }
#   cent
# }


kernel_vec<-function(x,y, kernel_fn, squared=FALSE,...){
  n<-nrow(as.matrix(x))
  m<-nrow(as.matrix(y))
  k_orig<-k_plain(x,y,kernel_fn = bm_kernel_L2)
  K_gram_matrix<-k_plain(x,y=NULL, bm_kernel_L2)
  first_comp<-as.matrix(colSums(K_gram_matrix))
  second_comp<-colSums(k_orig)
  third_comp<-sum(K_gram_matrix)
  #cent<-k_orig-1/n*first_comp-1/n*second_comp+1/(n^2)*third_comp
  cent<-apply(k_orig, 2, function(col) col - 1/n*first_comp)
  cent<-apply(cent,1, function(row) row - 1/n*second_comp)
  cent<-cent+1/(n^2)*third_comp
  cent<- t(matrix(cent,nrow = m))
  if (squared==TRUE) {
    cent=K_gram(x,y=NULL, bm_kernel_L2, centred=TRUE, squared=FALSE)%*%(cent)
  }
  cent
}



