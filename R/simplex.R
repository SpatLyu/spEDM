methods::setGeneric("simplex", function(data, ...) standardGeneric("simplex"))

.simplex_sf_method = \(data,target,lib,pred = lib,E = 1:10,tau = 1,k = E+2,
                       nb = NULL, threads = detectThreads(), trend.rm = TRUE){
  vec = .uni_lattice(data,target,trend.rm)
  lib = .check_indices(lib,length(vec))
  pred = .check_indices(pred,length(vec))
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  res = RcppSimplex4Lattice(vec,nb,lib,pred,E,k,tau,threads)
  return(.bind_xmapself(res))
}

.simplex_spatraster_method = \(data,target,lib,pred = lib,E = 1:10,tau = 1,k = E+2,
                               threads = detectThreads(), trend.rm = TRUE){
  mat = .uni_grid(data,target,trend.rm)
  res = RcppSimplex4Grid(mat,lib,pred,E,k,tau,threads)
  return(.bind_xmapself(res))
}

#' simplex forecast
#'
#' @inheritParams embedded
#' @param lib Row numbers(`vector` for lattice data) or row-column numbers(`matrix` for grid data) to create the library from observations.
#' @param pred (optional) Row numbers(`vector` for lattice data) or row-column numbers(`matrix` for grid data) used for predictions.
#' @param k (optional) Number of nearest neighbors used for prediction.
#' @param threads (optional) Number of threads.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{self mapping prediction results}
#' }
#' @export
#'
#' @name simplex
#' @rdname simplex
#' @aliases simplex,sf-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData"))
#' \donttest{
#' simplex(columbus,target = "CRIME",lib = 1:49)
#' }
methods::setMethod("simplex", "sf", .simplex_sf_method)

#' @rdname simplex
methods::setMethod("simplex", "SpatRaster", .simplex_spatraster_method)
