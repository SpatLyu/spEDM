methods::setGeneric("simplex", function(data, ...) standardGeneric("simplex"))

.simplex_sf_method = \(data,target,lib = NULL,pred = NULL,E = 1:10,tau = 1,k = E+2,
                       nb = NULL, threads = detectThreads(), detrend = TRUE){
  vec = .uni_lattice(data,target,detrend)
  if (is.null(lib)) lib = which(!is.na(vec))
  if (is.null(pred)) pred = lib
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  res = RcppSimplex4Lattice(vec,nb,lib,pred,E,k,tau,threads)
  return(.bind_xmapself(res,target))
}

.simplex_spatraster_method = \(data,target,lib = NULL,pred = NULL,E = 1:10,tau = 1,
                               k = E+2, threads = detectThreads(), detrend = TRUE){
  mat = .uni_grid(data,target,detrend)
  if (is.null(lib)) lib = which(!is.na(mat), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  res = RcppSimplex4Grid(mat,lib,pred,E,k,tau,threads)
  return(.bind_xmapself(res,target))
}

#' simplex forecast
#'
#' @inheritParams embedded
#' @param lib (optional) Libraries indices.
#' @param pred (optional) Predictions indices.
#' @param k (optional) Number of nearest neighbors used in prediction.
#' @param threads (optional) Number of threads.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{self mapping prediction results}
#' \item{\code{varname}}{name of target variable}
#' }
#' @export
#'
#' @name simplex
#' @rdname simplex
#' @aliases simplex,sf-method
#' @references
#' Sugihara G. and May R. 1990. Nonlinear forecasting as a way of distinguishing chaos from measurement error in time series. Nature, 344:734-741.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' simplex(columbus,"crime")
#' }
methods::setMethod("simplex", "sf", .simplex_sf_method)

#' @rdname simplex
methods::setMethod("simplex", "SpatRaster", .simplex_spatraster_method)
