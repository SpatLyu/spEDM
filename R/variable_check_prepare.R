.check_character = \(x,...){
  xstrs = c(x,...)
  for (i in xstrs){
    if (!inherits(i,"character")) {
      stop("Please check the characters in the function input.")
    }
  }
  return(xstrs)
}

.check_inputelementnum = \(x,n){
  if (length(x) == 1) x = rep(x,times = n)
  if (length(x) > 2) x = x[seq(1,n,by = 1)]
  return(x)
}

.check_indices = \(x,totalnum){
  return(x[(x<=totalnum)&(x>=1)])
}

.uni_lattice = \(data,target){
  target = .check_character(target)
  res = data[,target,drop = TRUE]
  return(res)
}

.uni_grid = \(data,target){
  target = .check_character(target)
  data = data[[target]]
  res = matrix(terra::values(data),nrow = terra::nrow(data),byrow = TRUE)
  return(res)
}
