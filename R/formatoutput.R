.internal_xmapdf_print = \(x,keyname = "libsizes",significant = FALSE){
  resdf = x[[1]]
  bidirectional = x[[3]]
  if (bidirectional){
    if (significant) {
      resdf = resdf[resdf$x_xmap_y_sig < 0.05 & resdf$y_xmap_x_sig < 0.05,
                    c(keyname, "y_xmap_x_mean", "x_xmap_y_mean")]
    } else {
      resdf = resdf[,c(keyname, "y_xmap_x_mean", "x_xmap_y_mean")]
    }
    names(resdf) = c(keyname,
                     paste0(x$varname[1], "->", x$varname[2]),
                     paste0(x$varname[2], "->", x$varname[1]))
  } else {
    if (significant) {
      resdf = resdf[resdf$y_xmap_x_sig < 0.05,c(keyname, "y_xmap_x_mean")]
    } else {
      resdf = resdf[,c(keyname, "y_xmap_x_mean")]
    }
    names(resdf) = c(keyname,
                     paste0(x$varname[1], "->", x$varname[2]))
  }
  return(resdf)
}

#' print ccm result
#' @noRd
#' @export
print.ccm_res = \(x,...){
  print(.internal_xmapdf_print(x))
}

#' print cmc result
#' @noRd
#' @export
print.cmc_res = \(x,...){
  print(.internal_xmapdf_print(x,"neighbors",TRUE))
}

#' print pcm result
#' @noRd
#' @export
print.pcm_res = \(x,...){
  pxmap = x[-2]
  xmap = x[-1]

  cat('-------------------------------------- \n')
  cat("***partial cross mapping prediction*** \n")
  cat('-------------------------------------- \n')
  print(.internal_xmapdf_print(pxmap))
  cat("\n------------------------------ \n")
  cat("***cross mapping prediction*** \n")
  cat('------------------------------ \n')
  print(.internal_xmapdf_print(xmap))
}

#' print xmap_self result
#' @noRd
#' @export
print.xmap_self = \(x,...){
  res = x$xmap
  if (ncol(res) == 5){
    outres = OptEmbedDim(res)
    cat(paste0("The suggested E and k for variable ",x$varname," is ",outres[1]," and ",outres[2]), "\n")
  } else if (ncol(res) == 4){
    cat(paste0("The suggested theta for variable ",x$varname," is ",OptThetaParm(res)), "\n")
  } else {
    print(x)
  }
}

#' print sc result
#' @noRd
#' @export
print.sc_res = \(x,...){
  sc = round(x$sc,3)
  varname = x$varname
  cat(c("spatial causality test",
        paste0(paste0(x$varname[1], " -> ", x$varname[2],": statistic = "), sc[1], ", p value = ", sc[2]),
        paste0(paste0(x$varname[2], " -> ", x$varname[1],": statistic = "), sc[3], ", p value = ", sc[4])),
      sep = "\n")
}

#' plot ccm result
#' @noRd
#' @export
plot.ccm_res = \(x, family = "serif",
                 legend_texts = NULL,
                 legend_cols = c("#608dbe","#ed795b"),
                 xbreaks = NULL, xlimits = NULL,
                 ybreaks = seq(0, 1, by = 0.1),
                 ylimits = c(-0.05, 1), ...){
  resdf = x[[1]]
  bidirectional = x[[3]]

  if(is.null(xbreaks)) xbreaks = resdf$libsizes
  if(is.null(xlimits)) xlimits = c(min(xbreaks)-1,max(xbreaks)+1)
  if (is.null(legend_texts)) legend_texts = c(paste0(x$varname[1], " xmap ", x$varname[2]),
                                              paste0(x$varname[2], " xmap ", x$varname[1]))
  legend_texts = .check_inputelementnum(legend_texts,2)
  legend_cols = .check_inputelementnum(legend_cols,2)
  names(legend_cols) = c("x xmap y","y xmap x")

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
                                expand = c(0, 0), name = "Library size") +
    ggplot2::scale_y_continuous(breaks = ybreaks, limits = ylimits,
                                expand = c(0, 0), name = expression(rho)) +
    ggplot2::scale_color_manual(values = legend_cols,
                                labels = legend_texts,
                                name = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(family = family),
                   axis.text.x = ggplot2::element_text(angle = 30),
                   axis.title = ggplot2::element_text(family = family),
                   panel.grid = ggplot2::element_blank(),
                   legend.position = "inside",
                   legend.justification = c(0.05,1),
                   legend.background = ggplot2::element_rect(fill = 'transparent'),
                   legend.text = ggplot2::element_text(family = family))
  return(fig1)
}

#' plot pcm result
#' @noRd
#' @export
plot.pcm_res = \(x, partial = TRUE, ...){
  indice = ifelse(partial,-2,-1)
  xmap = x[indice]
  class(xmap) = "ccm"
  fig1 = plot.ccm_res(xmap,...)
  return(fig1)
}
