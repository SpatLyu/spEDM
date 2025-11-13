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

.check_distmetric = \(dist.metric){
  dm = 2
  if (dist.metric != "L2"){
    dm = 1
  }
  return(dm)
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

.internal_grid2df = \(data, grid.coord = TRUE){
  if (grid.coord) {
    dtf = terra::as.data.frame(data,xy = TRUE,na.rm = FALSE)
  } else {
    dtf = terra::as.data.frame(data,xy = FALSE,na.rm = FALSE)
    cellnum = terra::rowColFromCell(data,seq_len(terra::ncell(data)))
    colnames(cellnum) = c("x","y")
    dtf = dplyr::bind_cols(dtf,cellnum)
  }
  return(dtf)
}

.internal_detrend = \(data,.varname,coords = NULL){
  if (!is.null(coords)) {
    data = dplyr::bind_cols(data, coords)
    formula_suffix = "~X+Y"
  } else {
    formula_suffix = "~x+y"
  }

  for (.var in .varname){
    data[[.var]] = sdsfun::rm_lineartrend(paste0(.var,formula_suffix), data = data)
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

.uni_grid = \(data,target,detrend = FALSE,grid.coord = TRUE){
  if (is.null(target)) return(matrix(0,terra::nrow(data),terra::ncol(data)))
  target = .check_character(target)
  data = data[[target]]
  names(data) = "target"
  dtf = .internal_grid2df(data,grid.coord)
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

.multivar_grid = \(data,columns,detrend = FALSE,grid.coord = TRUE){
  columns = .check_character(columns)
  data = data[[columns]]
  .varname = paste0("z",seq_along(columns))
  names(data) = .varname
  dtf = .internal_grid2df(data,grid.coord)
  if (detrend){
    dtf = .internal_detrend(dtf,.varname)
  }
  res = as.matrix(dtf[,.varname,drop = FALSE])
  return(res)
}

.internal_xmapdf_binding = \(x_xmap_y, y_xmap_x, bidirectional,
                             keyname = "libsizes", only_cs = FALSE){
  if (only_cs){
    colnames(y_xmap_x) = c(keyname,"y_xmap_x_mean")
  } else {
    colnames(y_xmap_x) = c(keyname,"y_xmap_x_mean","y_xmap_x_sig",
                           "y_xmap_x_lower","y_xmap_x_upper")
  }
  y_xmap_x = as.data.frame(y_xmap_x)

  if (bidirectional){
    if (only_cs){
      colnames(x_xmap_y) = c(keyname,"x_xmap_y_mean")
    } else {
      colnames(x_xmap_y) = c(keyname,"x_xmap_y_mean","x_xmap_y_sig",
                             "x_xmap_y_lower","x_xmap_y_upper")
    }
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
  resdf = .internal_xmapdf_binding(x_xmap_y,y_xmap_x,bidirectional)
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

  txmap = .internal_xmapdf_binding(txxmapy,tyxmapx,bidirectional)
  dxmap = .internal_xmapdf_binding(dxxmapy,dyxmapx,bidirectional)

  res = list("pxmap" = dxmap, "xmap" = txmap,
             "varname" = varname[1:2],
             "bidirectional" = bidirectional)
  class(res) = 'pcm_res'
  return(res)
}

.bind_intersectdf = \(varname,x_xmap_y,y_xmap_x,bidirectional){
  xmapdf = .internal_xmapdf_binding(x_xmap_y$xmap,y_xmap_x$xmap,bidirectional,keyname = "neighbors")
  csdf = .internal_xmapdf_binding(x_xmap_y$cs,y_xmap_x$cs,bidirectional,only_cs = TRUE)
  res = list("xmap" = xmapdf, "cs" = csdf, "varname" = varname, "bidirectional" = bidirectional)
  class(res) = 'cmc_res'
  return(res)
}

.bind_xmapself = \(x,varname,method,tau = NULL,...){
  res = list("xmap" = as.data.frame(x),"varname" = varname,"method" = method)
  if (!is.null(tau)) res = append(res,c("tau" = tau))
  class(res) = "xmap_self"
  return(res)
}

.bind_sc = \(sc,varname,...){
  res = list("sc" = sc,"varname" = varname)
  class(res) = "sc_res"
  return(res)
}

.bind_slm = \(mat_list,x,y,z,transient){
  if (is.null(transient)) {
    res = lapply(mat_list, \(.x) apply(.x,1,mean,na.rm = TRUE))
  } else {
    res = lapply(mat_list, \(.x) apply(.x[,-unique(abs(transient)),drop = FALSE],1,mean,na.rm = TRUE))
  }

  indices = NULL
  if (is.null(x)) indices = c(indices,1)
  if (is.null(y)) indices = c(indices,2)
  if (is.null(z)) indices = c(indices,3)
  if (is.null(indices)) return(res)
  return(res[-indices])
}

.run_gpc = \(x, y, E, k, tau, style, lib, pred, dist.metric, zero.tolerance, relative,
             weighted, na.rm, threads, bidirectional = FALSE, varname = NULL, nb = NULL,
             libsizes = NULL, boot = 9, random = TRUE, seed = 42, parallel.level = "low",
             progressbar = FALSE){
  dist.metric = .check_distmetric(dist.metric)
  if (is.null(libsizes)){
    if (is.null(nb)){
      res = RcppGPC4Grid(y,x,lib,pred,E,tau,style,k,zero.tolerance,dist.metric,relative,weighted,na.rm,threads)
    } else {
      res = RcppGPC4Lattice(y,x,nb,lib,pred,E,tau,style,k,zero.tolerance,dist.metric,relative,weighted,na.rm,threads)
    }
    res$causality$direction = "y_xmap_x"
    res$summary$direction = "y_xmap_x"
    res$pattern = list("y_xmap_x" = as.data.frame(res$pattern))
    res$varname = varname
    res$bidirectional = FALSE

    if (bidirectional){
      res_bi = .run_gpc(y, x, E, k, tau, style, lib, pred,
                        dist.metric, zero.tolerance,
                        relative, weighted, na.rm, threads,
                        FALSE, varname, nb)
      res_bi$causality$direction = "x_xmap_y"
      res_bi$summary$direction = "x_xmap_y"
      res$causality = rbind(res$causality,res_bi$causality)
      res$summary = rbind(res$summary,res_bi$summary)
      res$pattern$x_xmap_y = res_bi$pattern$y_xmap_x
      res$bidirectional = TRUE
    }

    class(res) = "pc_res"
    return(res)
  } else {
    pl = .check_parallellevel(parallel.level)
    if (is.null(nb)){
      res = RcppGPCRobust4Grid(y,x,libsizes,lib,pred,E,tau,style,k,boot,random,seed,zero.tolerance,
                               dist.metric, relative, weighted, na.rm, threads, pl, progressbar)
    } else {
      res = RcppGPCRobust4Lattice(y,x,nb,libsizes,lib,pred,E,tau,style,k,boot,random,seed,zero.tolerance,
                                  dist.metric, relative, weighted, na.rm, threads, pl, progressbar)
    }
    res$direction = "y_xmap_x"

    if (bidirectional){
      res_bi = .run_gpc(y, x, E, k, tau, style, lib, pred, dist.metric,
                        zero.tolerance, relative, weighted, na.rm, threads, FALSE,
                        varname, nb, libsizes, boot, random, seed, pl, progressbar)
      res_bi$xmap$direction = "x_xmap_y"
      res$xmap = rbind(res$xmap,res_bi$xmap)
    }

    return(structure(list(xmap = res, varname = varname,
                          bidirectional = bidirectional),
                     class = "rpc_res"))
  }
}
