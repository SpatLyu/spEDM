#' geographical convergent cross mapping
#'
#' @param cause Name of causal variable.
#' @param effect Name of effect variable.
#' @param data The observation data, must be `sf` or `SpatRaster` object.
#' @param libsizes A vector of library sizes to use.
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
#' \donttest{
#' gccm("HOVAL", "CRIME", data = columbus, libsizes = seq(5,45,5))
#' }
gccm = \(cause, effect, data, libsizes, E = 3,
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
    if (trendRM){
      # cause = RcppLinearTrendRM(cause,as.double(coords[,1]),as.double(coords[,2]))
      # effect = RcppLinearTrendRM(effect,as.double(coords[,1]),as.double(coords[,2]))
      dtf = data.frame(cause = cause, effect = effect, x = coords[,1], y = coords[,2])
      cause = sdsfun::rm_lineartrend("cause~x+y", data = dtf)
      effect = sdsfun::rm_lineartrend("effect~x+y", data = dtf)
    }

    y_causes_x = RcppGCCM4Lattice(cause,effect,nb,libsizes,E)
    x_causes_y = RcppGCCM4Lattice(effect,cause,nb,libsizes,E)

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
    selvec = seq(5,maxlibsize,5)
    if (is.null(RowCol)) RowCol = as.matrix(expand.grid(selvec,selvec))

    y_causes_x = RcppGCCM4Grid(causemat,effectmat,libsizes,RowCol,E)
    x_causes_y = RcppGCCM4Grid(effectmat,causemat,libsizes,RowCol,E)

  } else {
    stop("The data should be `sf` or `SpatRaster` object!")
  }

  colnames(y_causes_x) = c("lib_sizes","y_causes_x_mean","y_causes_x_sig",
                           "y_causes_x_upper","y_causes_x_lower")
  y_causes_x = as.data.frame(y_causes_x)
  colnames(x_causes_y) = c("lib_sizes","x_causes_y_mean","x_causes_y_sig",
                           "x_causes_y_upper","x_causes_y_lower")
  x_causes_y = as.data.frame(x_causes_y)

  res = y_causes_x |>
    dplyr::left_join(x_causes_y, by = "lib_sizes") |>
    dplyr::arrange(lib_sizes)
  return(res)
}
