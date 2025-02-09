methods::setGeneric("embedded", function(data, ...) standardGeneric("embedded"))

.embedded_sf_method = \(data,target,E = 3,tau = 1,nb = NULL){
  vec = .uni_lattice(data,target)
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  return(RcppGenLatticeEmbeddings(vec,nb,E,tau))
}

.embedded_spatraster_method = \(data,target,E = 3,tau = 1){
  mat = .uni_grid(data,target)
  return(RcppGenGridEmbeddings(mat,E,tau))
}

#' generate embeddings
#'
#' @param data The observation data.
#' @param target Name of target variable.
#' @param E (optional) The dimensions of the embedding.
#' @param tau (optional) The step of spatial lags.
#' @param nb (optional) The neighbours list.
#'
#' @return A matrix
#' @export
#'
#' @name embedded
#' @rdname embedded
#' @aliases embedded,sf-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData"))
#' embedded(columbus,target = "CRIME", E = 3)
#'
methods::setMethod("embedded", "sf", .embedded_sf_method)

#' @rdname embedded
methods::setMethod("embedded", "SpatRaster", .embedded_spatraster_method)
