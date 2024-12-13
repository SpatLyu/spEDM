#' geographical convergent cross mapping
#'
#' @param cause Name of causal variable.
#' @param effect Name of effect variable.
#' @param data The observation data, must be `sf` or `SpatRaster` object.
#' @param libsizes (optional) A vector of library sizes to use.
#' @param E (optional) The dimensions of the embedding.
#' @param nb (optional) The neighbours list.
#' @param RowCol (optional) Matrix of selected row and cols numbers.
#' @param trendRM (optional) Whether to remove the linear trend.
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData")[1],
#'                        quiet=TRUE)
#' gccm("HOVAL", "CRIME", data = columbus)
gccm = \(cause, effect, data, libsizes = NULL, E = 3,
         nb = NULL, RowCol = NULL, trendRM =TRUE) {
  if (!inherits(cause,"character") || !inherits(effect,"character")) {
    stop("The `cause` and `effect` must be character.")
  }

  if (inherits(data,"sf")) {
    coords = sdsfun::sf_coordinates(data)
    cause = data[,cause,drop = TRUE]
    effect = data[,effect,drop = TRUE]
    if (is.null(nb)) nb = sdsfun::spdep_nb(data)
    if (length(cause) != length(nb)) stop("Incompatible Data Dimensions!")
    if (is.null(libsizes)) libsizes = floor(seq(E + 2,length(cause),
                                                length.out = floor(sqrt(length(cause)))))
    if (trendRM){
      # cause = RcppLinearTrendRM(cause,as.double(coords[,1]),as.double(coords[,2]))
      # effect = RcppLinearTrendRM(effect,as.double(coords[,1]),as.double(coords[,2]))
      dtf = data.frame(cause = cause, effect = effect, x = coords[,1], y = coords[,2])
      cause = sdsfun::rm_lineartrend("cause~x+y", data = dtf)
      effect = sdsfun::rm_lineartrend("effect~x+y", data = dtf)
    }

    x_xmap_y = RcppGCCM4Lattice(cause,effect,nb,libsizes,E)
    y_xmap_x = RcppGCCM4Lattice(effect,cause,nb,libsizes,E)

  } else if (inherits(data,"SpatRaster")) {
    data = data[[c(cause,effect)]]
    names(data) = c("cause","effect")

    dtf = terra::as.data.frame(data,xy = TRUE,na.rm = FALSE)
    if (trendRM){
      dtf$cause = sdsfun::rm_lineartrend("cause~x+y", data = dtf)
      dtf$effect = sdsfun::rm_lineartrend("effect~x+y", data = dtf)
    }
    causemat = matrix(dtf[,3],nrow = terra::nrow(data),byrow = TRUE)
    effectmat = matrix(dtf[,4],nrow = terra::nrow(data),byrow = TRUE)

    maxlibsize = min(dim(causemat))
    if (is.null(libsizes)) libsizes = floor(seq(E + 2, maxlibsize,
                                                length.out = floor(sqrt(maxlibsize))))
    selvec = seq(5,maxlibsize,5)
    if (is.null(RowCol)) RowCol = as.matrix(expand.grid(selvec,selvec))

    x_xmap_y = RcppGCCM4Grid(causemat,effectmat,libsizes,RowCol,E)
    y_xmap_x = RcppGCCM4Grid(effectmat,causemat,libsizes,RowCol,E)

  } else {
    stop("The data should be `sf` or `SpatRaster` object!")
  }

  colnames(x_xmap_y) = c("lib_sizes","x_xmap_y_mean","x_xmap_y_sig",
                         "x_xmap_y_upper","x_xmap_y_lower")
  x_xmap_y = as.data.frame(x_xmap_y)
  colnames(y_xmap_x) = c("lib_sizes","y_xmap_x_mean","y_xmap_x_sig",
                         "y_xmap_x_upper","y_xmap_x_lower")
  y_xmap_x = as.data.frame(y_xmap_x)

  res = x_xmap_y |>
    dplyr::left_join(y_xmap_x, by = "lib_sizes") |>
    dplyr::arrange(lib_sizes)
  return(res)
}
