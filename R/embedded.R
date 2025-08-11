.embedded_sf_method = \(data,target,E = 3,tau = 1,style = 1,nb = NULL,detrend = FALSE){
  vec = .uni_lattice(data,target,detrend)
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  return(RcppGenLatticeEmbeddings(vec,nb,E,tau,style))
}

.embedded_spatraster_method = \(data,target,E = 3,tau = 1,style = 1,detrend = FALSE){
  mat = .uni_grid(data,target,detrend)
  return(RcppGenGridEmbeddings(mat,E,tau,style))
}

#' embedding spatial cross sectional data
#'
#' @param data observation data.
#' @param target name of target variable.
#' @param E (optional) embedding dimensions.
#' @param tau (optional) step of spatial lags.
#' @param style (optional) embedding style; 0 includes current state, 1 excludes it.
#' @param nb (optional) neighbours list.
#' @param detrend (optional) whether to remove the linear trend.
#'
#' @return A matrix
#' @export
#' @name embedded
#' @aliases embedded,sf-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' v = embedded(columbus,"crime")
#' v[1:5,]
#'
#' cu = terra::rast(system.file("case/cu.tif", package="spEDM"))
#' r = embedded(cu,"cu")
#' r[1:5,]
#'
methods::setMethod("embedded", "sf", .embedded_sf_method)

#' @rdname embedded
methods::setMethod("embedded", "SpatRaster", .embedded_spatraster_method)
