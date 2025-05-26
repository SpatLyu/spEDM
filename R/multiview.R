methods::setGeneric("multiview", function(data, ...) standardGeneric("multiview"))

.multiview_sf_method = \(data,column,target,nvar,lib = NULL,pred = NULL,E = 3,tau = 1,k = E+2,
                         nb = NULL, top = NULL, threads = detectThreads(), detrend = TRUE){
  xmat = .multivar_lattice(data,column,detrend)
  yvec = .uni_lattice(data,target,detrend)
  if (is.null(lib)) lib = .internal_library(cbind(xmat,yvec))
  if (is.null(pred)) pred = lib
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  if (is.null(top)) top = 0
  res = RcppMultiView4Lattice(xmat,yvec,nb,lib,pred,E,tau,k,top,nvar,threads)
  return(res)
}

.multiview_spatraster_method = \(data,column,target,nvar,lib = NULL,pred = NULL,E = 3,tau = 1,
                                 k = E+2,top = NULL,threads = detectThreads(),detrend = TRUE){
  xmat = .multivar_grid(data,column,detrend)
  ymat = .multivar_grid(data,target,detrend)
  if (is.null(lib)) lib = .internal_library(cbind(xmat,ymat),TRUE)
  if (is.null(pred)) pred = lib
  if (is.null(top)) top = 0
  res = RcppMultiView4Grid(xmat,ymat,lib,pred,E,tau,k,top,nvar,threads)
  return(res)
}

#' multiview embedding forecast
#'
#' @inheritParams simplex
#' @param nvar Number of variable combinations.
#' @param top (optional) Number of reconstructions used in MVE forecast.
#'
#' @return A vector (when input is sf object) or matrix
#' @export
#'
#' @name multiview
#' @rdname multiview
#' @aliases multiview,sf-method
#' @references
#' Ye H., and G. Sugihara, 2016. Information leverage in interconnected ecosystems: Overcoming the curse of dimensionality. Science 353:922-925.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' multiview(columbus,
#'           column = c("inc","crime","open","plumb","discbd"),
#'           target = "hoval", nvar = 3)
#' }
methods::setMethod("multiview", "sf", .multiview_sf_method)

#' @rdname multiview
methods::setMethod("multiview", "SpatRaster", .multiview_spatraster_method)
