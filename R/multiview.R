methods::setGeneric("multiview", function(data, ...) standardGeneric("multiview"))

.multiview_sf_method = \(data,columns,target,nvar,lib,pred = lib,E = 3,tau = 1,k = 4,
                         nb = NULL, top = NULL, threads = detectThreads(), trend.rm = TRUE){
  xmat = .multivar_lattice(data,columns,trend.rm)
  yvec = .uni_lattice(data,target,trend.rm)
  lib = .check_indices(lib,length(yvec))
  pred = .check_indices(pred,length(yvec))
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  if (is.null(top)) top = 0
  res = RcppMultiView4Lattice(xmat,yvec,nb,lib,pred,E,tau,k,top,nvar,threads)
  return(res)
}

.multiview_spatraster_method = \(data,columns,target,nvar,lib,pred = lib,E = 3,tau = 1,k = 4,
                                 top = NULL, threads = detectThreads(), trend.rm = TRUE){
  xmat = .multivar_grid(data,columns,trend.rm)
  ymat = .uni_grid(data,target,trend.rm)
  if (is.null(top)) top = 0
  res = RcppMultiView4Grid(xmat,ymat,lib,pred,E,tau,k,top,nvar,threads)
  return(res)
}

#' multiview embedding
#'
#' @inheritParams simplex
#' @param columns Individual variable names.
#' @param nvar Number of variable combinations.
#' @param top (optional) Number of reconstructions used for MVE forecast.
#'
#' @return A vector (when input is sf object) or matrix
#' @export
#'
#' @name multiview
#' @rdname multiview
#' @aliases multiview,sf-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData"))
#' \donttest{
#' multiview(columbus,
#'           columns = c("INC","CRIME","OPEN","PLUMB","DISCBD"),
#'           target = "HOVAL", nvar = 3,lib = 1:39, pred = 40:49)
#' }
methods::setMethod("multiview", "sf", .multiview_sf_method)

#' @rdname multiview
methods::setMethod("multiview", "SpatRaster", .multiview_spatraster_method)
