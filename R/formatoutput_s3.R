#' print ccm result
#' @noRd
#' @export
print.ccm_res = \(x,...){
  resdf = x$xmap
  bidirectional = x$bidirectional

  if (bidirectional){
    resdf = resdf[,c("libsizes", "y_xmap_x_mean", "x_xmap_y_mean")]
    names(resdf) = c('libsizes',
                     paste0(x$varname[1], "->", x$varname[2]),
                     paste0(x$varname[2], "->", x$varname[1]))
  } else {
    resdf = resdf[,c("libsizes", "y_xmap_x_mean")]
    names(resdf) = c('libsizes',
                     paste0(x$varname[1], "->", x$varname[2]))
  }
  print(resdf)
}

#' print pcm result
#' @noRd
#' @export
print.pcm_res = \(x,...){
  pxmap = x[-2]
  xmap = x[-1]
  class(pxmap) = "ccm"
  class(xmap) = "ccm"

  cat("partial cross mapping prediction: \n")
  print.ccm_res(pxamp)
  cat("cross mapping prediction: \n")
  print.ccm_res(xamp)
}

#' plot ccm result
#' @noRd
#' @export
plot.ccm_res = \(x, family = "serif", xbreaks = NULL, xlimits = NULL,
                 ybreaks = seq(0, 1, by = 0.1), ylimits = c(-0.05, 1), ...){
  resdf = x$xmap
  bidirectional = x$bidirectional
  if (bidirectional){
    resdf = resdf[,c("libsizes", "x_xmap_y_mean", "y_xmap_x_mean")]
  } else {
    resdf = resdf[,c("libsizes", "y_xmap_x_mean")]
  }
  if(is.null(xbreaks)) xbreaks = resdf$libsizes
  if(is.null(xlimits)) xlimits = c(min(xbreaks)-1,max(xbreaks)+1)

  if (bidirectional){
    resdf = resdf[,c("libsizes", "x_xmap_y_mean", "y_xmap_x_mean")]
  } else {
    resdf = resdf[,c("libsizes", "y_xmap_x_mean")]
  }

  fig1 = ggplot2::ggplot(data = resdf,
                         ggplot2::aes(x = libsizes)) +
    ggplot2::geom_line(ggplot2::aes(y = y_xmap_x_mean,
                                    color = "y xmap x"),
                       lwd = 1.25)

  if (bidirectional){
    fig1 = fig1 + ggplot2::geom_line(ggplot2::aes(y = x_xmap_y_mean,
                                                  color = "x xmap y"),
                                     lwd = 1.25)
  }

  fig1 = fig1 +
    ggplot2::scale_x_continuous(breaks = xbreaks, limits = xlimits,
                                expand = c(0, 0), name = "Lib of Sizes") +
    ggplot2::scale_y_continuous(breaks = ybreaks, limits = ylimits,
                                expand = c(0, 0), name = expression(rho)) +
    ggplot2::scale_color_manual(values = c("x xmap y" = "#608dbe",
                                           "y xmap x" = "#ed795b"),
                                labels = c(paste0(x$varname[2], " causes ", x$varname[1]),
                                           paste0(x$varname[1], " causes ", x$varname[2])),
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

#' plot pcm result
#' @noRd
#' @export
plot.pcm_res = \(x, ...){
  pxmap = x[-2]
  xmap = x[-1]
  class(pxmap) = "ccm"
  class(xmap) = "ccm"

  cat("plot partial cross mapping prediction: \n")
  print.ccm_res(pxamp,...)
  cat("plot cross mapping prediction: \n")
  print.ccm_res(xamp,...)
}
