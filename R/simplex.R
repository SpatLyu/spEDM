.simplex_sf_method = \(data,column,target,lib = NULL,pred = NULL,E = 1:10,tau = 1,
                       k = E+2, nb = NULL, threads = detectThreads(), detrend = TRUE){
  vx = .uni_lattice(data,column,detrend)
  vy = .uni_lattice(data,target,detrend)
  if (is.null(lib)) lib = .internal_library(cbind(vx,vy))
  if (is.null(pred)) pred = lib
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  res = RcppSimplex4Lattice(vx,vy,nb,lib,pred,E,k,tau,threads)
  return(.bind_xmapself(res,target))
}

.simplex_spatraster_method = \(data,column,target,lib = NULL,pred = NULL,E = 1:10,tau = 1,
                               k = E+2, threads = detectThreads(), detrend = TRUE){
  mx = .uni_grid(data,column,detrend)
  my = .uni_grid(data,target,detrend)
  if (is.null(lib)) lib = which(!(is.na(mx) | is.na(my)), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  res = RcppSimplex4Grid(mx,my,lib,pred,E,k,tau,threads)
  return(.bind_xmapself(res,target))
}

#' simplex forecast
#'
#' @inheritParams embedded
#' @param column name of library variable.
#' @param lib (optional) libraries indices.
#' @param pred (optional) predictions indices.
#' @param k (optional) number of nearest neighbors used in prediction.
#' @param threads (optional) number of threads to use.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{self mapping prediction results}
#' \item{\code{varname}}{name of target variable}
#' }
#' @export
#' @name simplex
#' @aliases simplex,sf-method
#' @references
#' Sugihara G. and May R. 1990. Nonlinear forecasting as a way of distinguishing chaos from measurement error in time series. Nature, 344:734-741.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' simplex(columbus,"inc","crime")
#' }
methods::setMethod("simplex", "sf", .simplex_sf_method)

#' @rdname simplex
methods::setMethod("simplex", "SpatRaster", .simplex_spatraster_method)
