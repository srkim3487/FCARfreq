#' @title Confidence Interval for the Spatial Dependence Parameter
#'
#' @description 
#' Computes the 95\% confidence interval for the spatial dependence parameter \eqn{\rho}.
#' 
#' @param p Selected truncation level.
#'
#' @param rhoHat Estimate of the spatial dependence parameter \eqn{\rho}.
#'
#' @param adj_mat An binary adjagency matrix defining the spatial neighborhood structure.
#' 
#' @return A numeric vector of length 2 giving the lower and upper bounds of the 95\% confidence interval for \eqn{\rho}.
#
#' @seealso
#' \code{\link{FCAR_est_fun}}
#' 
#' @examples
#' 
#' library(FCARfreq)
#' 
#' # Load example data included in the package
#' data(data_example)
#' data <- data_example$data
#' nbd_index <- data_example$nbd_index
#' 
#' # Create an adjacency matrix from nbd_index
#' n <- nrow(data)
#' adj_mat <- matrix(0, nrow = n, ncol = n)
#' for(i in 1:n){
#'  adj_mat[i, nbd_index[[i]]] <- 1
#'}
#'
#' # Estimation of alpha
#' alphaHat <- colMeans(data)
#' 
#'  est_FVE95 <- FCAR_est_fun(alphaHat, rhoHat=0.3, adj_mat, nbd_index, data, sel_p="FVE95")
#'  
#' CI_fun(est_FVE95$p, est_FVE95$rhoHat, adj_mat)
#' 
#' @export
CI_fun <- function(p, rhoHat, adj_mat){
  n <- dim(adj_mat)[1]
  # m <- dim(data)[2]
  # t <- seq(from = 0, to = 1, length.out = m)
  
  ## number of neigbors
  num_nei <- rowSums(adj_mat)
  
  cov_temp <- solve(diag(num_nei) - rhoHat *adj_mat) 
  l_sec_der_temp <- adj_mat %*% cov_temp
  neg_l_sec_der <- p / (2*n) * sum(diag(l_sec_der_temp %*% l_sec_der_temp))
  
  sd <- sqrt(n)*sqrt(neg_l_sec_der)
  CI <- c(rhoHat-qnorm(1 - 0.05/2)/sd, rhoHat+qnorm(1 - 0.05/2)/sd)
  return(CI)
}




