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

#' spatial logistic map
#'
#' @param data observation data.
#' @param x first spatial variable.
#' @param y (optional) second spatial variable.
#' @param z (optional) third spatial variable.
#' @param k (optional) number of neighbors to used.
#' @param step (optional) number of simulation time steps.
#' @param alpha_x (optional) growth parameter for x.
#' @param alpha_y (optional) growth parameter for y.
#' @param alpha_z (optional) growth parameter for y.
#' @param beta_xy (optional) cross-inhibition from x to y.
#' @param beta_xz (optional) cross-inhibition from x to z.
#' @param beta_yx (optional) cross-inhibition from y to x.
#' @param beta_yz (optional) cross-inhibition from y to z.
#' @param beta_zx (optional) cross-inhibition from z to x.
#' @param beta_zy (optional) cross-inhibition from z to y.
#' @param threshold (optional) set to `NaN` if the absolute value exceeds this threshold.
#'
#' @return A list
#' @export
#'
#' @name slm
#' @rdname slm
#' @aliases slm,sf-method
#' @references
#' Willeboordse, F.H., The spatial logistic map as a simple prototype for spatiotemporal chaos, Chaos, 533–540 (2003).
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' columbus$inc = sdsfun::normalize_vector(columbus$inc)
#' slm(columbus,"inc")
#'
methods::setMethod("slm", "sf", .slm_sf_method)

#' @rdname slm
methods::setMethod("slm", "SpatRaster", .slm_spatraster_method)
