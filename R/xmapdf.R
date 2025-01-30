.bind_xmapdf = \(varname,x_xmap_y,y_xmap_x,bidirectional = TRUE){

  colnames(y_xmap_x) = c("libsizes","y_xmap_x_mean","y_xmap_x_sig",
                         "y_xmap_x_upper","y_xmap_x_lower")
  y_xmap_x = as.data.frame(y_xmap_x)

  if (bidirectional){
    colnames(x_xmap_y) = c("libsizes","x_xmap_y_mean","x_xmap_y_sig",
                           "x_xmap_y_upper","x_xmap_y_lower")
    x_xmap_y = as.data.frame(x_xmap_y)
    resdf = x_xmap_y |>
      dplyr::full_join(y_xmap_x, by = "libsizes") |>
      dplyr::arrange(libsizes)
  } else {
    resdf = dplyr::arrange(y_xmap_x,libsizes)
  }

  res = list("xmap" = resdf, "varname" = varname, "bidirectional" = bidirectional)
  class(res) = 'ccm_res'
  return(res)
}
