methods::setGeneric("smap", function(data, ...) standardGeneric("smap"))

.smap_sf_method = \(data, target, lib, pred, E = 3, k = 4,
                    theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03,
                              0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8),
                    nb = NULL,threads = detectThreads(),include.self = FALSE){
  vec = .uni_lattice(data,target)
  lib = .check_indices(lib,length(vec))
  pred = .check_indices(pred,length(vec))
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  res = RcppSMap4Lattice(vec,nb,lib,pred,theta,E,k,threads,include.self)
  cat(paste0("The suggested theta for variable ",target," is ",OptThetaParm(res)), "\n")
  return(res)
}

.smap_spatraster_method = \(data, target, lib, pred, E = 3, k = 4,
                            theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03,
                                      0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8),
                            threads = detectThreads(), include.self = FALSE){
  mat = .uni_grid(data,target)
  res = RcppSMap4Grid(mat,lib,pred,theta,E,k,threads,include.self)
  cat(paste0("The suggested theta for variable ",target," is ",OptThetaParm(res)), "\n")
  return(res)
}

#' smap forecasting
#'
#' @inheritParams simplex
#' @theta (optional) Weighting parameter for distances
#'
#' @return A matrix
#' @export
#'
#' @name smap
#' @rdname smap
#' @aliases smap,sf-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData")[1],
#'                        quiet=TRUE)
#' \donttest{
#' smap(columbus,target = "CRIME",lib = 1:29,pred = 30:49)
#' }
methods::setMethod("smap", "sf", .smap_sf_method)

#' @rdname smap
methods::setMethod("smap", "SpatRaster", .smap_spatraster_method)
