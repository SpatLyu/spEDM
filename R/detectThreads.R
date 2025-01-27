#' detect the number of available threads
#'
#' @return An integer
#' @export
#'
#' @examples
#' \donttest{
#' detectThreads()
#' }
detectThreads = \() {
  return(DetectMaxNumThreads())
}
