#' Save plot as PDF
#' 
#' @description
#' `save_pdf()` function saves a plot as PDF.
#' 
#' @param p plot.
#' @inheritParams grDevices::pdf
#' @inheritDotParams grDevices::pdf
#' 
#' @seealso [grDevices::pdf()]
#' 
#' @export

save_pdf <- function(p, file, width, height, ...) {
  grDevices::pdf(file = file, width = width, height = height, ...)
  grDevices::pdf.options(encoding = 'CP1250')
  print(p)
  invisible(grDevices::dev.off())
}