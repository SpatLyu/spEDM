gccm = \(cause, effect, nbmat = NULL, data = NULL,
         libsizes = NULL, E = 3, ...) {
  if (is.null(nbmat)){
    if(inherits(data,"sf")){
      nbmat = spdep::nb2mat(sdsfun::spdep_nb(data,...),
                            style = "B", zero.policy = TRUE)
    } else {
      stop("When `nbmat` is NULL, the data must be provided as `sf` object!")
    }
  }

  if (inherits(cause,"character") & inherits(effect,"character")){
    if (is.null(data)){
      stop("When `cause` and `effect` are character, the data must be provided!")
    }
    cause = data[,cause,drop = TRUE]
    effect = data[,effect,drop = TRUE]
  }

  n = length(cause)
  if (n != nrow(nbmat)) stop("Incompatible Data Dimensions!")
  if (is.null(libsizes)) libsizes = floor(seq(1,n,length.out = 15))

  x_xmap_y = RcppGCCMLattice(cause,effect,nbmat,libsizes,E)
  return(x_xmap_y)
}
