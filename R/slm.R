methods::setGeneric("slm", function(data, ...) standardGeneric("slm"))

.slm_sf_method = \(data, x, y = NULL, z = NULL,
                   k = 4, step = 20, alpha_x = 0.77, alpha_y = 0, alpha_z = 0,
                   beta_xy = 0, beta_xz = 0, beta_yx = 0, beta_yz = 0, beta_zx = 0, beta_zy = 0,
                   threshold = 1e10, nb = NULL){
  vx = .uni_lattice(data,x,FALSE)
  vy = .uni_lattice(data,y,FALSE)
  vz = .uni_lattice(data,z,FALSE)
  if (is.null(nb)) nb = .internal_lattice_nb(data)
  return(lapply(RcppSLMTri4Lattice(vx,vy,vz,nb,k,step,alpha_x,alpha_y,alpha_z,beta_xy,beta_xz,beta_yx,beta_yz,beta_zx,beta_zy,threshold),
                \(.x) apply(.x,1,mean,na.rm = TRUE)))
}

.slm_spatraster_method = \(data, x, y = NULL, z = NULL,
                           k = 4, step = 20, alpha_x = 0.77, alpha_y = 0, alpha_z = 0,
                           beta_xy = 0, beta_xz = 0, beta_yx = 0, beta_yz = 0, beta_zx = 0, beta_zy = 0,
                           threshold = 1e10){
  mx = .uni_grid(data,x,FALSE)
  my = .uni_grid(data,y,FALSE)
  mz = .uni_grid(data,z,FALSE)
  return(lapply(RcppSLMTri4Grid(mx,my,mz,k,step,alpha_x,alpha_y,alpha_z,beta_xy,beta_xz,beta_yx,beta_yz,beta_zx,beta_zy,threshold),
                \(.x) apply(.x,1,mean,na.rm = TRUE)))
}
