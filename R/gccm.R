#' geographical convergent cross mapping
#'
#' @param cause Name of causal variable.
#' @param effect Name of effect variable.
#' @param data The observation data, must be `sf` or `SpatRaster` object.
#' @param libsizes A vector of library sizes to use.
#' @param E (optional) The dimensions of the embedding.
#' @param tau (optional) The step of spatial lags.
#' @param k (optional) Number of nearest neighbors to use for prediction.
#' @param nb (optional) The neighbours list.
#' @param RowCol (optional) Matrix of selected row and cols numbers.
#' @param trendRM (optional) Whether to remove the linear trend.
#' @param progressbar (optional) whether to print the progress bar.
#'
#' @return A list.
#' \describe{
#' \item{\code{xmap}}{cross-mapping prediction outputs}
#' \item{\code{varname}}{names of causal and effect variable}
#' }
#' @export
#'
#' @examples
#' columbus = sf::read_sf(system.file("shapes/columbus.gpkg", package="spData")[1],
#'                        quiet=TRUE)
#' \donttest{
#' g = gccm("HOVAL", "CRIME", data = columbus, libsizes = seq(5,45,5))
#' g
#' plot(g,ylimits = c(0,0.65))
#' }
gccm = \(cause, effect, data, libsizes, E = 3, tau = 1, k = E + 2,
         nb = NULL, RowCol = NULL, trendRM = TRUE, progressbar = TRUE) {
  if (!inherits(cause,"character") || !inherits(effect,"character")) {
    stop("The `cause` and `effect` must be character.")
  }
  varname = c(cause,effect)

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

    x_xmap_y = RcppGCCM4Lattice(cause,effect,nb,libsizes,E,tau,k,progressbar)
    y_xmap_x = RcppGCCM4Lattice(effect,cause,nb,libsizes,E,tau,k,progressbar)

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

    x_xmap_y = RcppGCCM4Grid(causemat,effectmat,libsizes,RowCol,E,tau,k,progressbar)
    y_xmap_x = RcppGCCM4Grid(effectmat,causemat,libsizes,RowCol,E,tau,k,progressbar)

  } else {
    stop("The data should be `sf` or `SpatRaster` object!")
  }

  colnames(x_xmap_y) = c("libsizes","x_xmap_y_mean","x_xmap_y_sig",
                         "x_xmap_y_upper","x_xmap_y_lower")
  x_xmap_y = as.data.frame(x_xmap_y)
  colnames(y_xmap_x) = c("libsizes","y_xmap_x_mean","y_xmap_x_sig",
                         "y_xmap_x_upper","y_xmap_x_lower")
  y_xmap_x = as.data.frame(y_xmap_x)

  resdf = x_xmap_y |>
    dplyr::left_join(y_xmap_x, by = "libsizes") |>
    dplyr::arrange(libsizes)

  res = list("xmap" = resdf, "varname" = varname)
  class(res) = 'ccm_res'
  return(res)
}
