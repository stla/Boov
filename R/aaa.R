#' @useDynLib Boov, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Objects imported from other packages
#' @description These objects are imported from other packages.
#'   Follow the links to their documentation:
#'   \code{\link[PolygonSoup:toRGL]{toRGL}},
#'   \code{\link[PolygonSoup:plotEdges]{plotEdges}}.
#' @importFrom PolygonSoup toRGL plotEdges
#' @export toRGL plotEdges
#' @name Boov-imports
#' @aliases toRGL plotEdges
#' @docType import
NULL

utils::globalVariables("exterior")
