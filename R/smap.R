methods::setGeneric("smap", function(data, ...) standardGeneric("smap"))

.smap_sf_method = \(data, target, lib, pred = lib, E = 3, tau = 1, k = 4,
                    theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03,
                              0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8),
                    nb = NULL, threads = detectThreads(), trend.rm = TRUE){
  vec = .uni_lattice(data,target,trend.rm)
  lib = .check_indices(lib,length(vec))
  pred = .check_indices(pred,length(vec))
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  res = RcppSMap4Lattice(vec,nb,lib,pred,theta,E,tau,k,threads)
  cat(paste0("The suggested theta for variable ",target," is ",OptThetaParm(res)), "\n")
  return(res)
}

.smap_spatraster_method = \(data, target, lib, pred = lib, E = 3, tau = 1, k = 4,
                            theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03,
                                      0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8),
                            threads = detectThreads(), trend.rm = TRUE){
  mat = .uni_grid(data,target,trend.rm)
  res = RcppSMap4Grid(mat,lib,pred,theta,E,tau,k,threads)
  cat(paste0("The suggested theta for variable ",target," is ",OptThetaParm(res)), "\n")
  return(res)
}

#' smap forecasting
#'
#' @inheritParams simplex
#' @param theta (optional) Weighting parameter for distances
#'
#' @return A matrix
#' @export
#'
#' @name smap
#' @rdname smap
#' @aliases smap,sf-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData"))
#' \donttest{
#' smap(columbus,target = "INC",lib = 1:49)
#' }
methods::setMethod("smap", "sf", .smap_sf_method)

#' @rdname smap
methods::setMethod("smap", "SpatRaster", .smap_spatraster_method)
