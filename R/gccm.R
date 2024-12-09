#' geographical convergent cross mapping
#'
#' @param cause The causal series (as a character or a vector).
#' @param effect The effect series (as a character or a vector).
#' @param nbmat (optional) The neighbours matrix.
#' @param data (optional) The observation data.
#' @param libsizes (optional) A vector of library sizes to use
#' @param E (optional) The number of dimensions for the attractor reconstruction.
#' @param ... (optional) Other parameters passed to `sdsfun::spdep_nb()`.
#'
#' @return A `tibble`.
#' @export
#'
#' @examples
#' columbus <- sf::read_sf(system.file("shapes/columbus.gpkg", package="spData")[1],
#'                         quiet=TRUE)
#' gccm("CRIME","HOVAL", data = columbus)
gccm = \(cause, effect, nbmat = NULL, data = NULL,
         libsizes = NULL, E = 3, ...) {
  if (is.null(nbmat)){
    if(inherits(data,"sf")){
      nbmat = spdep::nb2mat(sdsfun::spdep_nb(data,...),
                            style = "B", zero.policy = TRUE)
    } else {
      stop("When `nbmat` is NULL, the data must be provided as `sf` object!")
    }
  }

  if (inherits(cause,"character") & inherits(effect,"character")){
    if (is.null(data)){
      stop("When `cause` and `effect` are character, the data must be provided!")
    }
    cause = data[,cause,drop = TRUE]
    effect = data[,effect,drop = TRUE]
  }

  if (length(cause) != nrow(nbmat)) stop("Incompatible Data Dimensions!")
  if (is.null(libsizes)) libsizes = floor(seq(1,length(cause),length.out = 15))

  x_xmap_y = RcppGCCMLattice(cause,effect,nbmat,libsizes,E)
  y_xmap_x = RcppGCCMLattice(effect,cause,nbmat,libsizes,E)
  return(x_xmap_y)
}
