methods::setGeneric("embedded", function(data, ...) standardGeneric("embedded"))

.embedded_sf_method = \(data,target,E = 3,tau = 1,nb = NULL,trend.rm = FALSE){
  vec = .uni_lattice(data,target,trend.rm)
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  return(RcppGenLatticeEmbeddings(vec,nb,E,tau))
}

.embedded_spatraster_method = \(data,target,E = 3,tau = 1,trend.rm = FALSE){
  mat = .uni_grid(data,target,trend.rm)
  return(RcppGenGridEmbeddings(mat,E,tau))
}

#' embedding spatial cross sectional data
#'
#' @param data The observation data.
#' @param target Name of target variable.
#' @param E (optional) Dimensions of the embedding.
#' @param tau (optional) Step of spatial lags.
#' @param nb (optional) The neighbours list.
#' @param trend.rm (optional) Whether to remove the linear trend.
#'
#' @return A matrix
#' @export
#'
#' @name embedded
#' @rdname embedded
#' @aliases embedded,sf-method
#' @references
#'  \insertRef{Takens}{spEDM}
#'  \insertRef{Mane}{spEDM}
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' embedded(columbus,target = "crime", E = 3)
#'
methods::setMethod("embedded", "sf", .embedded_sf_method)

#' @rdname embedded
methods::setMethod("embedded", "SpatRaster", .embedded_spatraster_method)
