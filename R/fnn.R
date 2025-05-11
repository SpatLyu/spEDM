methods::setGeneric("fnn", function(data, ...) standardGeneric("fnn"))

.fnn_sf_method = \(data, target, lib = NULL, pred = NULL, E = 1:10, tau = 1, nb = NULL,
                   rt = 10, eps = 2, threads = detectThreads(), trend.rm = TRUE){
  vec = .uni_lattice(data,target,trend.rm)
  rt = .check_inputelementnum(rt,length(E))
  eps = .check_inputelementnum(eps,length(E))
  if (is.null(lib)) lib = which(!is.na(vec))
  if (is.null(pred)) pred = lib
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  return(RcppFNN4Lattice(vec,nb,rt,eps,lib,pred,E,tau,threads))
}

.fnn_spatraster_method = \(data, target, lib = NULL, pred = NULL, E = 1:10, tau = 1,
                           rt = 10, eps = 2, threads = detectThreads(), trend.rm = TRUE){
  mat = .uni_grid(data,target,trend.rm)
  rt = .check_inputelementnum(rt,length(E))
  eps = .check_inputelementnum(eps,length(E))
  if (is.null(lib)) lib = which(!is.na(mat), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  return(RcppFNN4Grid(mat,rt,eps,lib,pred,E,tau,threads))
}

#' false nearest neighbours
#'
#' @inheritParams embedded
#' @param lib (optional) Libraries indices.
#' @param pred (optional) Predictions indices.
#' @param rt (optional) escape factor.
#' @param eps (optional) neighborhood diameter.
#' @param threads (optional) Number of threads.
#'
#' @return A vector
#' @export
#'
#' @name fnn
#' @rdname fnn
#' @aliases fnn,sf-method
#' @references
#' Kennel M. B., Brown R. and Abarbanel H. D. I., Determining embedding dimension for phase-space reconstruction using a geometrical construction, Phys. Rev. A, Volume 45, 3403 (1992).
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' fnn(columbus,"crime")
#' }
methods::setMethod("fnn", "sf", .fnn_sf_method)

#' @rdname fnn
methods::setMethod("fnn", "SpatRaster", .fnn_spatraster_method)

