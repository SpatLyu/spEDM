methods::setGeneric("multiview", function(data, ...) standardGeneric("multiview"))

.multiview_sf_method = \(data,columns,target,nvar,lib = NULL,pred = NULL,E = 3,tau = 1,k = E+2,
                         nb = NULL, top = NULL, threads = detectThreads(), trend.rm = TRUE){
  xmat = .multivar_lattice(data,columns,trend.rm)
  yvec = .uni_lattice(data,target,trend.rm)
  if (is.null(lib)) lib = 1:nrow(data)
  if (is.null(pred)) pred = lib
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  if (is.null(top)) top = 0
  res = RcppMultiView4Lattice(xmat,yvec,nb,lib,pred,E,tau,k,top,nvar,threads)
  return(res)
}

.multiview_spatraster_method = \(data,columns,target,nvar,lib = NULL,pred = NULL,E = 3,tau = 1,
                                 k = E+2,top = NULL,threads = detectThreads(),trend.rm = TRUE){
  xmat = .multivar_grid(data,columns,trend.rm)
  ymat = .uni_grid(data,target,trend.rm)
  if (is.null(lib)) lib = .internal_samplemat(ymat)
  if (is.null(pred)) pred = lib
  if (is.null(top)) top = 0
  res = RcppMultiView4Grid(xmat,ymat,lib,pred,E,tau,k,top,nvar,threads)
  return(res)
}

#' multiview embedding forecast
#'
#' @inheritParams simplex
#' @param columns Names of individual variables.
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
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' multiview(columbus,
#'           columns = c("inc","crime","open","plumb","discbd"),
#'           target = "hoval", nvar = 3)
#' }
methods::setMethod("multiview", "sf", .multiview_sf_method)

#' @rdname multiview
methods::setMethod("multiview", "SpatRaster", .multiview_spatraster_method)
