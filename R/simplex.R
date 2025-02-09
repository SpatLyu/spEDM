methods::setGeneric("simplex", function(data, ...) standardGeneric("simplex"))

.simplex_sf_method = \(data,target,lib,pred = lib,E = 1:10,k = 4,nb = NULL,
                       threads = detectThreads(),include.self = FALSE){
  vec = .uni_lattice(data,target)
  lib = .check_indices(lib,length(vec))
  pred = .check_indices(pred,length(vec))
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  res = RcppSimplex4Lattice(vec,nb,lib,pred,E,k,threads,include.self)
  cat(paste0("The suggested embedding dimension E for variable ",target," is ",OptEmdedDim(res)), "\n")
  return(res)
}

.simplex_spatraster_method = \(data,target,lib,pred = lib,E = 1:10,k = 4,
                               threads = detectThreads(),include.self = FALSE){
  mat = .uni_grid(data,target)
  res = RcppSimplex4Grid(mat,lib,pred,E,k,threads,include.self)
  cat(paste0("The suggested embedding dimension E for variable ",target," is ",OptEmdedDim(res)), "\n")
  return(res)
}

#' simplex forecasting
#'
#' @param data The observation data.
#' @param target Name of target variable.
#' @param lib The row numbers(`vector`) of lattice data or the row-column numbers(`matrix`) of grid data for creating the library from observations.
#' @param pred (optional) The row numbers(`vector`) of lattice data or the row-column numbers(`matrix`) of grid data used for predictions.
#' @param E (optional) The dimensions of the embedding.
#' @param k (optional) Number of nearest neighbors to use for prediction.
#' @param nb (optional) The neighbours list.
#' @param threads (optional) Number of threads.
#' @param include.self (optional) Whether to include the current state when constructing the embedding vector.
#'
#' @return A matrix
#' @export
#'
#' @name simplex
#' @rdname simplex
#' @aliases simplex,sf-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData"))
#' \donttest{
#' simplex(columbus,target = "CRIME",lib = 1:29,pred = 30:49)
#' }
methods::setMethod("simplex", "sf", .simplex_sf_method)

#' @rdname simplex
methods::setMethod("simplex", "SpatRaster", .simplex_spatraster_method)
