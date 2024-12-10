#' geographical convergent cross mapping
#'
#' @param cause The causal series (as a character or a vector).
#' @param effect The effect series (as a character or a vector).
#' @param nb (optional) The neighbours list.
#' @param coords (optional) The coordinates matrix.
#' @param data (optional) The observation data.
#' @param libsizes (optional) A vector of library sizes to use.
#' @param E (optional) The dimensions of the embedding.
#' @param ... (optional) Other parameters passed to `sdsfun::spdep_nb()`.
#'
#' @return A `tibble`.
#' @export
#'
#' @examples
#' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData")[1],
#'                        quiet=TRUE)
#' gccm("HOVAL", "CRIME", data = columbus)
gccm = \(cause, effect, nb = NULL, coords = NULL,
         data = NULL, libsizes = NULL, E = 3, ...) {

  if (is.null(nb)){
    if(inherits(data,"sf")){
      nb = sdsfun::spdep_nb(data,...)
    } else {
      stop("When `nb` is NULL, the data must be provided as an `sf` object!")
    }
  }

  if (is.null(coords)){
    if(inherits(data,"sf")){
      coords = sdsfun::sf_coordinates(data)
    } else {
      stop("When `coords` is NULL, the data must be provided as an `sf` object!")
    }
  }
  if (inherits(coords,"character")) coords = as.matrix(data[,coords])

  if (inherits(cause,"character") || inherits(effect,"character")){
    if (is.null(data)){
      stop("When `cause` and `effect` are character, the data must be provided!")
    }
    cause = data[,cause,drop = TRUE]
    effect = data[,effect,drop = TRUE]
  }

  if (length(cause) != length(nb)) stop("Incompatible Data Dimensions!")
  if (is.null(libsizes)) libsizes = floor(seq(E + 2,length(cause),
                                              length.out = floor(sqrt(length(cause)))))

  # effect = RcppLinearTrendRM(effect,as.double(coords[,1]),as.double(coords[,2]))
  # cause = RcppLinearTrendRM(cause,as.double(coords[,1]),as.double(coords[,2]))
  dtf = data.frame(cause = cause, effect = effect,
                   x = coords[,1], y = coords[,2])
  effect = sdsfun::rm_lineartrend("effect~x+y", data = dtf)
  cause = sdsfun::rm_lineartrend("cause~x+y", data = dtf)

  x_xmap_y = RcppGCCMLattice(cause,effect,nb,libsizes,E)
  colnames(x_xmap_y) = c("lib_sizes","x_xmap_y_mean","x_xmap_y_sig",
                         "x_xmap_y_upper","x_xmap_y_lower")
  x_xmap_y = tibble::as_tibble(x_xmap_y)

  y_xmap_x = RcppGCCMLattice(effect,cause,nb,libsizes,E)
  colnames(y_xmap_x) = c("lib_sizes","y_xmap_x_mean","y_xmap_x_sig",
                         "y_xmap_x_upper","y_xmap_x_lower")
  y_xmap_x = tibble::as_tibble(y_xmap_x)

  res = x_xmap_y |>
    dplyr::left_join(y_xmap_x, by = "lib_sizes") |>
    dplyr::arrange(lib_sizes)
  return(res)
}
