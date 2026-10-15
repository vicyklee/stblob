# Generate threeblobs.rda -----------------------------------------------------

## function: gen_gaussian_data ------------------------------------------------
# n: an integer of the number of clusters
# size: an integer vector of sizes of a cluster. If an integer value is supplied, all clusters will be of this size.
# center: a data matrix or data frame of coordinates.
# ...:  the optional argument includes `sigma`.
# sigma: a covariance matrix or a list of covariance matrices to specify different covariance matrices. See [mvtnorm::rmvnorm()].

gen_gaussian_data <- function(n, size, center = NULL, ...) {
  args <- list(...)
  sigma <- args$sigma
  
  if (length(size) == 1) {
    size <- rep(size, n)
  }
  
  if (is.null(center)) {
    center <- matrix(sample(x = seq(-10,10,0.01),
                            size = 2*n,
                            replace = T),
                     nrow = n)
  }
  
  data <- list()
  for (i in 1:n) {
    if (is.null(sigma)) {
      data[[i]] <- mvtnorm::rmvnorm(n = size[i], mean = center[i,])
    } else {
      if (!is.list(sigma)) {
        data[[i]] <- mvtnorm::rmvnorm(n = size[i],
                                      mean = center[i,],
                                      sigma = sigma)
      } else {
        data[[i]] <- mvtnorm::rmvnorm(n = size[i],
                                      mean = center[i,],
                                      sigma = sigma[[i]])
      }
    }
    # label cluster
    data[[i]] <- cbind(data[[i]],i)
  }
  
  data <- do.call(rbind, data)
  colnames(data) <- c("x","y","clust")
  data <- as.data.frame(data)
  data$clust <- as.factor(data$clust)
  return(data)
}

## generate toy data ----------------------------------------------------------
# install.packages("mvtnorm")

set.seed(12)
x_y_type <- gen_gaussian_data(
  n = 3,
  size = c(100, 100, 100),
  center = as.matrix(data.frame(x = c(-10,-5,10),
                                y = c(-3,3,0))),
  sigma = matrix(c(2,0,0,2), nrow = 2)
)

# use "clust" as "type"
names(x_y_type)[3] <- "type"

z <- c(seq(1,300,1)[-seq(1,300,3)][1:100],
       seq(1,300,1)[-seq(1,300,3)][101:200],
       seq(1,300,3))

threeblobs <- cbind(sample_id = 1:nrow(x_y_type), x_y_type, z)
threeblobs <- as.data.frame(threeblobs[c("sample_id", "x", "y", "z", "type")])

usethis::use_data(threeblobs, overwrite = TRUE)
