#' detect the number of available threads
#'
#' @return An integer
#' @export
#'
#' @examples
#' spEDM::detectThreads()
#' 
detectThreads = \() {
  return(DetectMaxNumThreads())
}
