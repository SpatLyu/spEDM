#' geographical convergent cross mapping
#'
#' @param cause The causal series (as a character or a vector).
#' @param effect The effect series (as a character or a vector).
#' @param nbmat (optional) The neighbours matrix.
#' @param coords (optional) The coordinates matrix.
#' @param data (optional) The observation data.
#' @param libsizes (optional) A vector of library sizes to use.
#' @param E (optional) The dimensions of the embedding.
#' @param trend_rm (optional) Variables that need to remove linear trends.
#' @param ... (optional) Other parameters passed to `sdsfun::spdep_nb()`.
#'
#' @return A `tibble`.
#' @export
#'
#' @examples
#' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData")[1],
#'                        quiet=TRUE)
#' gccm("CRIME","HOVAL", data = columbus)
gccm = \(cause, effect, nbmat = NULL, coords = NULL, data = NULL, libsizes = NULL,
         E = 3, trend_rm = c("effect", "cause", "none", "all"), ...) {
  if (is.null(nbmat)){
    if(inherits(data,"sf")){
      nbmat = spdep::nb2mat(sdsfun::spdep_nb(data,...),
                            style = "B", zero.policy = TRUE)
    } else {
      stop("When `nbmat` is NULL, the data must be provided as `sf` object!")
    }
  }

  if (is.null(coords)){
    if(inherits(data,"sf")){
      coords = sdsfun::sf_coordinates(data)
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
  if (is.null(libsizes)) libsizes = floor(seq(E + 2,length(cause),length.out = 15))

  trend_rm = match.arg(trend_rm)
  switch(trend_rm,
         "effect" = {
           effect = RcppLinearTrendRM(effect,as.double(coords[,1]),as.double(coords[,2]))
         },
         "cause" = {
           cause = RcppLinearTrendRM(cause,as.double(coords[,1]),as.double(coords[,2]))
         },
         "all" = {
           effect = RcppLinearTrendRM(effect,as.double(coords[,1]),as.double(coords[,2]))
           cause = RcppLinearTrendRM(cause,as.double(coords[,1]),as.double(coords[,2]))
         })


  x_xmap_y = RcppGCCMLattice(cause,effect,nbmat,libsizes,E)
  colnames(x_xmap_y) = c("lib_sizes","x_xmap_y_mean","x_xmap_y_sig",
                         "x_xmap_y_upper","x_xmap_y_lower")
  x_xmap_y = tibble::as_tibble(x_xmap_y)

  y_xmap_x = RcppGCCMLattice(effect,cause,nbmat,libsizes,E)
  colnames(y_xmap_x) = c("lib_sizes","y_xmap_x_mean","y_xmap_x_sig",
                         "y_xmap_x_upper","y_xmap_x_lower")
  y_xmap_x = tibble::as_tibble(y_xmap_x)

  res = x_xmap_y |>
    dplyr::left_join(y_xmap_x, by = "lib_sizes") |>
    dplyr::arrange(lib_sizes)
  return(res)
}
