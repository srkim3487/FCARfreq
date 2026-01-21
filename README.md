The package implements the estimation and inference procedures proposed in:

Kim, S., (2026+). A new class of functional conditional autoregressive models.

To install and load the package:
```r
devtools::install_github("srkim3487/FCARfreq")
library(FCARfreq)
```

Below is an example using the included simulated dataset. For details on each function, please refer to the corresponding help page.

```r
# Load example data included in the package
data(data_example)
data <- data_example$data
nbd_index <- data_example$nbd_index

# Create an adjacency matrix from nbd_index
n <- nrow(data)
adj_mat <- matrix(0, nrow = n, ncol = n)
for(i in 1:n){
 adj_mat[i, nbd_index[[i]]] <- 1
}

# Estimation of alpha
alphaHat <- colMeans(data)

# Estimation
est_FVE95 <- FCAR_est_fun(alphaHat, rhoHat=0.3, adj_mat, nbd_index, data, sel_p="FVE95")

# Confidence Interval for the Spatial Dependence Parameter
CI_fun(est_FVE95$p, est_FVE95$rhoHat, adj_mat)
```
