methods::setGeneric("embed", function(data, ...) standardGeneric("embed"))

.embed_sf_method = \(data,target,E = 3,nb = NULL,include.self = FALSE){
  target = .check_tgcharacter(target)
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  scsvec = data[,target,drop = TRUE]
  return(RcppGenLatticeEmbeddings(scsvec,nb,E,include.self))
}

.embed_spatraster_method = \(data,target,E = 3,include.self = TRUE){
  target = .check_tgcharacter(target)
  data = data[[target]]
  mat = matrix(terra::values(data),nrow = terra::nrow(data),byrow = TRUE)
  return(RcppGenGridEmbeddings(mat,E,include.self))
}

#' generate embeddings
#'
#' @param data The observation data.
#' @param target Name of target variable.
#' @param E (optional) The dimensions of the embedding.
#' @param nb (optional) The neighbours list.
#' @param include.self (optional) Whether to include the current state.
#'
#' @return A matrix
#' @export
#'
#' @examples
#' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData")[1],
#'                        quiet=TRUE)
#' embed(columbus,target = "CRIME", E = 3)
#'
methods::setMethod("embed", "sf", .embed_sf_method)

#' @rdname embed
methods::setMethod("embed", "SpatRaster", .embed_spatraster_method)
