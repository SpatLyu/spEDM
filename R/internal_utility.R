.check_character = \(x,...){
  xstrs = c(x,...)
  for (i in xstrs){
    if (!inherits(i,"character")) {
      stop("Please check the characters in the function input.")
    }
  }
  return(xstrs)
}

.check_inputelementnum = \(x,n,condsnum = NULL){
  if (is.null(condsnum) || length(x) == 1){
    res = rep(x,length.out = n)
  } else if (length(x) == 2) {
    res = c(rep(x[1],2),rep(x[2],condsnum))
  } else {
    res = c(x[1:2],rep(x[c(-1,-2)],length.out = condsnum))
  }
  return(res)
}

.check_parallellevel = \(parallel.level){
  pl = 0
  if (parallel.level != "low"){
    pl = 1
  }
  return(pl)
}

.internal_varname = \(conds = NULL){
  .varname = c("cause","effect")
  if (!is.null(conds)){
    .varname = c(.varname,paste0("z",seq_along(conds)))
  }
  return(.varname)
}

.internal_lattice_nb = \(data){
  if (sdsfun::sf_geometry_type(data) %in% c('point','multipoint')){
    nb = sdsfun::spdep_nb(data,k = 8)
  } else {
    nb = sdsfun::spdep_nb(data)
  }
  return(nb)
}

.internal_detrend = \(data,.varname,coords = NULL){
  if (is.null(coords)){
    for (i in seq_along(.varname)){
      data[,.varname[i]] = sdsfun::rm_lineartrend(paste0(.varname[i],"~x+y"), data = data)
    }
  } else {
    data = dplyr::bind_cols(data,coords)
    for (i in seq_along(.varname)){
      data[,.varname[i]] = sdsfun::rm_lineartrend(paste0(.varname[i],"~X+Y"), data = data)
    }
  }
  return(data)
}

.internal_library = \(data,mat = FALSE){
  nnaindice = which(apply(is.na(data),1,\(.x) !any(.x)))
  if (mat) nnaindice = matrix(nnaindice,ncol = 1)
  return(nnaindice)
}

.uni_lattice = \(data,target,detrend = FALSE){
  if (is.null(target)) return(rep(0,nrow(data)))
  target = .check_character(target)
  coords = as.data.frame(sdsfun::sf_coordinates(data))
  data = sf::st_drop_geometry(data)
  data = data[,target,drop = FALSE]
  names(data) = "target"
  if (detrend){
    data = .internal_detrend(data,"target",coords)
  }
  res = data[,"target",drop = TRUE]
  return(res)
}

.uni_grid = \(data,target,detrend = FALSE){
  if (is.null(target)) return(matrix(0,terra::nrow(data),terra::ncol(data)))
  target = .check_character(target)
  data = data[[target]]
  names(data) = "target"
  dtf = terra::as.data.frame(data,xy = TRUE,na.rm = FALSE)
  if (detrend){
    dtf = .internal_detrend(dtf,"target")
  }
  res = matrix(dtf[,"target"],nrow = terra::nrow(data),byrow = TRUE)
  return(res)
}

.multivar_lattice = \(data,columns,detrend = FALSE){
  columns = .check_character(columns)
  coords = as.data.frame(sdsfun::sf_coordinates(data))
  data = sf::st_drop_geometry(data)
  data = data[,columns,drop = FALSE]
  .varname = paste0("z",seq_along(columns))
  names(data) = .varname
  if (detrend){
    data = .internal_detrend(data,.varname,coords)
  }
  res = as.matrix(data[,.varname,drop = FALSE])
  return(res)
}

.multivar_grid = \(data,columns,detrend = FALSE){
  columns = .check_character(columns)
  data = data[[columns]]
  .varname = paste0("z",seq_along(columns))
  names(data) = .varname
  dtf = terra::as.data.frame(data,xy = TRUE,na.rm = FALSE)
  if (detrend){
    dtf = .internal_detrend(dtf,.varname)
  }
  res = as.matrix(dtf[,.varname,drop = FALSE])
  return(res)
}

.internal_xmapdf_binding = \(varname,x_xmap_y,y_xmap_x,bidirectional,keyname = "libsizes"){
  colnames(y_xmap_x) = c(keyname,"y_xmap_x_mean","y_xmap_x_sig",
                         "y_xmap_x_upper","y_xmap_x_lower")
  y_xmap_x = as.data.frame(y_xmap_x)

  if (bidirectional){
    colnames(x_xmap_y) = c(keyname,"x_xmap_y_mean","x_xmap_y_sig",
                           "x_xmap_y_upper","x_xmap_y_lower")
    x_xmap_y = as.data.frame(x_xmap_y)
    resdf = x_xmap_y |>
      dplyr::full_join(y_xmap_x, by = keyname) |>
      dplyr::arrange({{keyname}})
  } else {
    resdf = dplyr::arrange(y_xmap_x,{{keyname}})
  }

  return(resdf)
}

.bind_xmapdf = \(varname,x_xmap_y,y_xmap_x,bidirectional){
  resdf = .internal_xmapdf_binding(varname,x_xmap_y,y_xmap_x,bidirectional)
  res = list("xmap" = resdf, "varname" = varname, "bidirectional" = bidirectional)
  class(res) = 'ccm_res'
  return(res)
}

.bind_xmapdf2 = \(varname,x_xmap_y,y_xmap_x,bidirectional){

  tyxmapx = y_xmap_x[,c(1,2,4:6),drop = FALSE]
  dyxmapx = y_xmap_x[,c(1,3,7:9),drop = FALSE]
  txxmapy = NULL
  dxxmapy = NULL
  if(bidirectional){
    txxmapy = x_xmap_y[,c(1,2,4:6),drop = FALSE]
    dxxmapy = x_xmap_y[,c(1,3,7:9),drop = FALSE]
  }

  txmap = .internal_xmapdf_binding(varname[1:2],txxmapy,tyxmapx,bidirectional)
  dxmap = .internal_xmapdf_binding(varname[1:2],dxxmapy,dyxmapx,bidirectional)

  res = list("pxmap" = dxmap, "xmap" = txmap,
             "varname" = varname[1:2],
             "bidirectional" = bidirectional)
  class(res) = 'pcm_res'
  return(res)
}

.bind_intersectdf = \(varname,x_xmap_y,y_xmap_x,bidirectional){
  resdf = .internal_xmapdf_binding(varname,x_xmap_y,y_xmap_x,bidirectional,keyname = "neighbors")
  res = list("xmap" = resdf, "varname" = varname, "bidirectional" = bidirectional)
  class(res) = 'cmc_res'
  return(res)
}

.bind_xmapself = \(x,varname,...){
  res = list("xmap" = x,"varname" = varname)
  class(res) = "xmap_self"
  return(res)
}

.bind_sc = \(sc,varname,...){
  res = list("sc" = sc,"varname" = varname)
  class(res) = "sc_res"
  return(res)
}

.bind_slm = \(mat_list,x,y,z,transient){
  res = lapply(mat_list, \(.x) apply(.x[,-transient,drop = FALSE],1,mean,na.rm = TRUE))
  indices = NULL
  if (is.null(x)) indices = c(indices,1)
  if (is.null(y)) indices = c(indices,2)
  if (is.null(z)) indices = c(indices,3)
  return(res[-indices])
}
