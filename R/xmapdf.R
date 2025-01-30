.bind_xmapdf = \(varname,y_xmap_x,x_xmap_y = NULL){
  bidirectional = TRUE
  colnames(x_xmap_y) = c("libsizes","x_xmap_y_mean","x_xmap_y_sig",
                         "x_xmap_y_upper","x_xmap_y_lower")
  x_xmap_y = as.data.frame(x_xmap_y)

  if (is.null(y_xmap_x)){
    resdf = dplyr::arrange(x_xmap_y,libsizes)
    bidirectional = FALSE
  } else {
    colnames(y_xmap_x) = c("libsizes","y_xmap_x_mean","y_xmap_x_sig",
                           "y_xmap_x_upper","y_xmap_x_lower")
    y_xmap_x = as.data.frame(y_xmap_x)

    resdf = x_xmap_y |>
      dplyr::left_join(y_xmap_x, by = "libsizes") |>
      dplyr::arrange(libsizes)
  }

  res = list("xmap" = resdf, "varname" = varname, "bidirectional" = TRUE)
  class(res) = 'ccm_res'
  return(res)
}
