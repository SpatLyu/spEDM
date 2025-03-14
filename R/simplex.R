methods::setGeneric("simplex", function(data, ...) standardGeneric("simplex"))

.simplex_sf_method = \(data,target,lib = NULL,pred = NULL,E = 1:10,tau = 1,k = E+2,
                       nb = NULL, threads = detectThreads(), trend.rm = TRUE){
  vec = .uni_lattice(data,target,trend.rm)
  if (is.null(lib)) lib = 1:nrow(data)
  if (is.null(pred)) pred = lib
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  res = RcppSimplex4Lattice(vec,nb,lib,pred,E,k,tau,threads)
  return(.bind_xmapself(res,target))
}

.simplex_spatraster_method = \(data,target,lib = NULL,pred = NULL,E = 1:10,tau = 1,
                               k = E+2, threads = detectThreads(), trend.rm = TRUE){
  mat = .uni_grid(data,target,trend.rm)
  if (is.null(lib)) lib = .internal_samplemat(mat)
  if (is.null(pred)) pred = lib
  res = RcppSimplex4Grid(mat,lib,pred,E,k,tau,threads)
  return(.bind_xmapself(res,target))
}

#' simplex forecast
#'
#' @inheritParams embedded
#' @param lib (optional) Libraries indices.
#' @param pred (optional) Predictions indices.
#' @param k (optional) Number of nearest neighbors used for prediction.
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
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' simplex(columbus,target = "crime")
#' }
methods::setMethod("simplex", "sf", .simplex_sf_method)

#' @rdname simplex
methods::setMethod("simplex", "SpatRaster", .simplex_spatraster_method)
