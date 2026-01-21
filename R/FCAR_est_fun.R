#' @title Estimation of Functional Conditional Autoregressive (FCAR) Model
#'
#' @description 
#' Implements the profile-based estimation procedure described in Algorithm 1 of the main paper.
#'
#' @param alphaHat Estimated value of \eqn{\alpha}, the mean function of the joint distribution.
#'
#' @param rhoHat An initial value of \eqn{\rho}, the spatial dependence parameter.
#'
#' @param adj_mat An binary adjagency matrix defining the spatial neighborhood structure.
#'
#' @param nbd_index An nb class object specifying the neighborhood structures, in which each element lists the adjacent spatial locations.
#'
#' @param data An n by T matrix of functional data. Each row represents one functional observation evaluated on a common grid.
#' 
#' @param sel_p A character string specifying the threshold level for the fraction of variance explained (FVE) used to select the truncation level. 
#' Must be one of \code{"FVE75"}, \code{"FVE80"}, \code{"FVE85"}, \code{"FVE90"}, \code{FVE95}, or \code{"FVE99"}. 
#' 
#' @param rho_tol Convergence tolerance for the spatial dependence parameter \eqn{\rho}.
#' 
#' @param G_tol Convergence tolerance for the covariance operator.
#' 
#' @return A list with the following elements:
#' \describe{
#'  \item{\code{rhoHat}:}{Estimate of the spatial dependence parameter \eqn{\rho}.}
#'  \item{\code{GHat}:}{Estimate of covariance operator.}
#'  \item{\code{iter}:}{Number of iterations used in the profile-based estimation procedure.}
#'  \item{\code{p}:}{Selected truncation level.}
#' }
#' @seealso
#' \code{\link{CI_fun}}
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
#' # Example: estimated spatial dependence parameter
#' est_FVE95$rhoHat
#' 
#' @export
FCAR_est_fun <- function(alphaHat, rhoHat, adj_mat, nbd_index, data, 
sel_p = c("FVE75", "FVE80", "FVE85", "FVE90", "FVE95", "FVE99"),
rho_tol = 1e-8, G_tol = 1e-6){
  iter <- 0
  GHat <- 0
  repeat{
    iter <- iter + 1
    
    # if(sel_p %in% c("FVE75", "FVE80", "FVE85", "FVE90", "FVE95", "FVE99")){
    ## Estimation of G
    covHat_list <- covHat_fun(rhoHat, alphaHat, adj_mat, nbd_index, data, sel_p = sel_p)
    GHat_update <- covHat_list$GHat
    p <- covHat_list$p
    lambdaHat <- covHat_list$lambdaHat[1:p]
    phiHat <- covHat_list$phiHat[,(1:p)]
    
    ############ optimization: rho
    pre <- precompute(adj_mat, data, phiHat)
    log_op_eigen <- function(rho){-logli_rho_eigen(rho, pre, lambdaHat)}
    opt <- optimize(log_op_eigen, c(0, 1))
    rho_update <- opt$minimum
    
    
    if((rhoHat-rho_update)^2 < rho_tol & mean((GHat - GHat_update)^2) < G_tol){break}
    if(iter > 50){
      return(NULL)  # discard and generate new dataset
    }
    rhoHat <- rho_update
    GHat <- GHat_update
  }
  return(list(rhoHat=rhoHat, GHat=GHat, 
              iter=iter, p=p
              ))
}


logli_rho_eigen <- function(rho, pre, lambdaHat){
    n <- length(pre$Dvec)
    
    # logdet part
    logdetQ <- sum(log(pre$Dvec)) + sum(log1p(-rho * pre$lambdas))
    
    # quadratic form for all p
    quad <- pre$q0 - rho * pre$q1
    
    -1/(2*n) * sum( quad/lambdaHat + n*log(lambdaHat) - logdetQ )
}


precompute <- function(adj_mat, data, phiHat){
  n <- nrow(data)
  Dvec <- rowSums(adj_mat)
  
  t <- seq(0,1,length.out=ncol(data))
  yHat <- data %*% phiHat * diff(range(t)) / ncol(data)
  
  # (1) Eigenvalues of standardized adjacency
  Dinv2 <- 1/sqrt(Dvec)
  Astd <- (Dinv2 * adj_mat) * Dinv2
  eig <- eigen(Astd, only.values=TRUE)
  lambdas <- eig$values
  
  # (2) Precompute quadratic form components
  # q0_j = y^T D y
  q0 <- colSums( (Dvec * yHat) * yHat )
  
  # q1_j = y^T W y
  Wy <- adj_mat %*% yHat
  q1 <- colSums(yHat * Wy)
  
  list(yHat=yHat, Dvec=Dvec, lambdas=lambdas, q0=q0, q1=q1)
}





covHat_fun <- function(rhoHat, alphaHat, adj_mat, nbd_index, data,
                       sel_p = c("FVE75", "FVE80", "FVE85", "FVE90", "FVE95", "FVE99")){
  # adj_mat: adjacency matrix,
  # nbd_index: nbd list
  # data: n * m matrix
  n <- dim(data)[1]
  m <- dim(data)[2]
  t <- seq(from = 0, to = 1, length.out = m)
  
  
  ## number of neigbors
  num_nei <- rowSums(adj_mat)
  
  ## centered data
  ZHat <- center_fun(rhoHat, alphaHat, adj_mat, nbd_index, data)
  
  ## sample covariance
  Ghat <- cov(ZHat)
  
  ## Mercer's Theorem
  eig <- eigen(Ghat)
  lambda <- eig$values * diff(range(t)) / m
  phi <- eig$vectors / sqrt(diff(range(t) / m))
  
  FVE <- cumsum(lambda / sum(lambda))
  if(sel_p == "FVE75"){
    p <- which(FVE > 0.75)[1]
  }
  if(sel_p == "FVE80"){
    p <- which(FVE > 0.80)[1]
  }
  if(sel_p == "FVE85"){
    p <- which(FVE > 0.85)[1]
  }
  if(sel_p == "FVE90"){
    p <- which(FVE > 0.90)[1]
  }
  if(sel_p == "FVE95"){
    p <- which(FVE > 0.95)[1]
  }
  if(sel_p == "FVE99"){
    p <- which(FVE > 0.99)[1]
  }
  
  
  return(list(GHat = Ghat, lambdaHat=lambda, phiHat=phi, p=p))
}


center_fun <- function(rhoHat, alphaHat, adj_mat, nbd_index, data){
  n <- dim(data)[1]
  m <- dim(data)[2]
  
  ## number of neigbors
  num_nei <- rowSums(adj_mat)
  
  ## centered data
  centered_y <- matrix(0, nrow = n ,ncol = m)
  for(k in 1:n){
    y <- data[k,]
    nbd_id <- nbd_index[[k]]
    
    # mean function
    if(length(nbd_id) > 1){
      mu_fun <- alphaHat + rhoHat/num_nei[k] * colSums(data[nbd_id,]-alphaHat) 
    } 
    if(length(nbd_id) == 1){
      mu_fun <- alphaHat + rhoHat/num_nei[k] * sum(data[nbd_id,]-alphaHat) 
    }
    
    
    # centered data
    centered_y[k,] <- (y - mu_fun) * sqrt(num_nei[k])
  }
  return(centered_y)
}



