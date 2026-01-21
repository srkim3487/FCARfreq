#' Example Dataset: data_example
#' 
#' This is a simulated dataset on 30*30 lattice.
#' 
#' @format A list with 2 elements:
#' \describe{
#' \item{\code{data}:}{A simulated spatial functional dataset consisting of 900 curves observed on 30*30 lattice. 
#' Each curve is evaluated at 50 equally spaced points over the interval [0, 1].}
#' \item{\code{nbd_index}:}{An \code{nb} class object specifying the neighborhood structures, 
#' in which each element lists the adjacent spatial locations.}
#' }
#'
#' @usage data(data_example)
#' @keywords datasets
"data_example"