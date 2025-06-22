.ic_sf_method = \(data, column, target, lib = NULL, pred = NULL, E = 2:10, tau = 1, k = E+2,
                  nb = NULL, threads = detectThreads(), parallel.level = "low", detrend = FALSE){
  vx = .uni_lattice(data,column,detrend)
  vy = .uni_lattice(data,target,detrend)
  if (is.null(lib)) lib = .internal_library(cbind(vx,vy))
  if (is.null(pred)) pred = lib
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  pl = .check_parallellevel(parallel.level)
  res = RcppIC4Lattice(vx,vy,nb,lib,pred,E,k,tau,0,threads,pl)
  return(.bind_xmapself(res,target,"ic",tau))
}

.ic_spatraster_method = \(data, column, target, lib = NULL, pred = NULL, E = 2:10, tau = 1, k = E+2,
                          threads = detectThreads(), parallel.level = "low", detrend = FALSE){
  mx = .uni_grid(data,column,detrend)
  my = .uni_grid(data,target,detrend)
  if (is.null(lib)) lib = which(!(is.na(mx) | is.na(my)), arr.ind = TRUE)
  if (is.null(pred)) pred = lib
  pl = .check_parallellevel(parallel.level)
  res = RcppIC4Grid(mx,my,lib,pred,E,k,tau,0,threads,pl)
  return(.bind_xmapself(res,target,"ic",tau))
}

#' intersection cardinality
#'
#' @inheritParams simplex
#' @param parallel.level (optional) level of parallelism, `low` or `high`.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{cross mapping performance}
#' \item{\code{varname}}{name of target variable}
#' \item{\code{method}}{method of cross mapping}
#' \item{\code{tau}}{step of time lag}
#' }
#' @export
#' @name ic
#' @aliases ic,sf-method
#' @references
#' Tao, P., Wang, Q., Shi, J., Hao, X., Liu, X., Min, B., Zhang, Y., Li, C., Cui, H., Chen, L., 2023. Detecting dynamical causality by intersection cardinal concavity. Fundamental Research.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' ic(columbus,"hoval","crime", k = 25)
#' }
methods::setMethod("ic", "sf", .ic_sf_method)

#' @rdname ic
methods::setMethod("ic", "SpatRaster", .ic_spatraster_method)
