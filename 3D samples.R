source("Kernels.R")
library(mvtnorm)
library(RColorBrewer)
library(plgp)

levelpersp <- function(x, y, z, colors=topo.colors, ...) {
  ## getting the value of the midpoint
  zz <- (z[-1,-1] + z[-1,-ncol(z)] + z[-nrow(z),-1] + z[-nrow(z),-ncol(z)])/4
  ## calculating the breaks
  breaks <- hist(zz, plot=FALSE)$breaks
  ## cutting up zz
  cols <- colors(length(breaks)-1)
  zzz <- cut(zz, breaks=breaks, labels=cols)
  ## plotting
  persp(x, y, z, col=as.character(zzz), ...)
  ## return breaks and colors for the legend
  list(breaks=breaks, colors=cols)
}
# Define the range and number of points for x1 and x2
n_points <- 20
x1 <- seq(0, 3, length.out = n_points)
x2 <- seq(0, 3, length.out = n_points)

# Sample from the bivariate normal distribution
N <- 3  # Number of samples
Y <- draw_samples_3D(x1, x2, N, seed = 125, kernel_fn = bm_kernel, centred = TRUE, squared = TRUE)

# Create the 3D plots with colors based on z-values using persp
par(mfrow = c(1, 3))
levelpersp(x1, x2, (matrix(Y[1, ], ncol = n_points)), theta = -30, phi = 30, colors = topo.colors)#, #col = matrix(colors, ncol = n_points),
#main = "Centered ANOVA kernel")

levelpersp(x1, x2, (matrix(Y[2, ], ncol = n_points)), theta = -30, phi = 30, colors = topo.colors)
levelpersp(x1, x2, (matrix(Y[3, ], ncol = n_points)), theta = -30, phi = 30, colors = topo.colors)