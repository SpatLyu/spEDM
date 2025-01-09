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
#' gccm("HOVAL", "CRIME", data = columbus, libsizes = seq(5,45,5))
#' }
gccm = \(cause, effect, data, libsizes, E = 3,
         nb = NULL, RowCol = NULL, trendRM =TRUE) {
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

  resdf = x_xmap_y |>
    dplyr::left_join(y_xmap_x, by = "lib_sizes") |>
    dplyr::arrange(lib_sizes)

  res = list("xmap" = resdf, "varname" = varname)
  class(res) = 'gccm_res'
  return(res)
}

#' print gccm result
#' @noRd
#' @export
print.gccm_res = \(x,...){
  resdf = x$xmap
  resdf = resdf[,c("lib_sizes", "x_xmap_y_mean", "y_xmap_x_mean")]
  names(resdf) = c('libsizes',
                   paste0(x$varname[2], "==>", x$varname[1]),
                   paste0(x$varname[1], "==>", x$varname[2]))
  print(resdf)
}

#' plot gccm result
#' @noRd
#' @export
plot.gccm_res = \(x,breaks = NULL,limits = NULL,family = "serif",...){
  resdf = x$xmap
  resdf = resdf[,c("lib_sizes", "x_xmap_y_mean", "y_xmap_x_mean")]

  if(is.null(breaks)) breaks = resdf$lib_sizes
  if(is.null(limits)) limits = c(min(breaks)-3,max(breaks)+3)
  fig1 = ggplot2::ggplot(data = resdf,
                         ggplot2::aes(x = lib_sizes)) +
    ggplot2::geom_line(ggplot2::aes(y = x_xmap_y_mean,
                                    color = "x xmap y"),
                       lwd = 1.25) +
    ggplot2::geom_line(ggplot2::aes(y = y_xmap_x_mean,
                                    color = "y xmap x"),
                       lwd = 1.25) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                                limits = c(-0.05, 1), expand = c(0, 0),
                                name = expression(rho)) +
    ggplot2::scale_x_continuous(name = "Lib of Sizes",
                                breaks = breaks,
                                limits = limits, expand = c(0, 0)) +
    ggplot2::scale_color_manual(values = c("x xmap y" = "#608dbe",
                                           "y xmap x" = "#ed795b"),
                                labels = c(paste0(x$varname[2], "causes", x$varname[1]),
                                           paste0(x$varname[1], "causes", x$varname[2])),
                                name = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(family = family),
                   axis.text.x = ggplot2::element_text(angle = 30),
                   axis.title = ggplot2::element_text(family = family),
                   panel.grid = ggplot2::element_blank(),
                   legend.position = "inside",
                   legend.justification = c('right','top'),
                   legend.background = ggplot2::element_rect(fill = 'transparent'),
                   legend.text = ggplot2::element_text(family = family))
  return(fig1)
}
