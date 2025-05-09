methods::setGeneric("smap", function(data, ...) standardGeneric("smap"))

.smap_sf_method = \(data,target,lib = NULL,pred = NULL,E = 3,tau = 1,k = E+2,
                    theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03,
                              0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8),
                    nb = NULL, threads = detectThreads(), trend.rm = TRUE){
  vec = .uni_lattice(data,target,trend.rm)
  if (is.null(lib)) lib = .internal_library(data)
  if (is.null(pred)) pred = lib
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  res = RcppSMap4Lattice(vec,nb,lib,pred,theta,E,tau,k,threads)
  return(.bind_xmapself(res,target))
}

.smap_spatraster_method = \(data,target,lib = NULL,pred = NULL,E = 3,tau = 1,k = E+2,
                            theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03,
                                      0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8),
                            threads = detectThreads(), trend.rm = TRUE){
  mat = .uni_grid(data,target,trend.rm)
  if (is.null(lib)) lib = .internal_samplemat(mat)
  if (is.null(pred)) pred = lib
  res = RcppSMap4Grid(mat,lib,pred,theta,E,tau,k,threads)
  return(.bind_xmapself(res,target))
}

#' smap forecast
#'
#' @inheritParams simplex
#' @param theta (optional) Weighting parameter for distances.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{self mapping prediction results}
#' \item{\code{varname}}{name of target variable}
#' }
#' @export
#'
#' @name smap
#' @rdname smap
#' @aliases smap,sf-method
#' @references
#' Sugihara G. 1994. Nonlinear forecasting for the classification of natural time series. Philosophical Transactions: Physical Sciences and Engineering, 348 (1688):477-495.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' smap(columbus,"inc")
#' }
methods::setMethod("smap", "sf", .smap_sf_method)

#' @rdname smap
methods::setMethod("smap", "SpatRaster", .smap_spatraster_method)
