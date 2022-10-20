#' @title Compound of five tetrahedra
#'
#' @description Five tetrahedra in a pretty configuration. Each tetrahedron
#'   is centered at the origin.
#'
#' @format A list with three fields: the field \code{meshes} is a list of five
#'   elements, each one representing a tetrahedron by a list with two elements,
#'   the vertices and the faces; the field \code{rglmeshes} is the list of
#'   the five corresponding \strong{rgl} meshes; the field \code{gmpmeshes} is
#'   the same as \code{meshes} except that the vertices are \strong{gmp}
#'   rational numbers.
#'
"tetrahedraCompound"
