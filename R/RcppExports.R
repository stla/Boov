# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

Intersection_EK <- function(rmeshes, clean, normals, triangulate) {
    .Call(`_Boov_Intersection_EK`, rmeshes, clean, normals, triangulate)
}

Intersection_Q <- function(rmeshes, clean, normals, triangulate) {
    .Call(`_Boov_Intersection_Q`, rmeshes, clean, normals, triangulate)
}

Difference_EK <- function(rmesh1, rmesh2, clean, normals, triangulate1, triangulate2) {
    .Call(`_Boov_Difference_EK`, rmesh1, rmesh2, clean, normals, triangulate1, triangulate2)
}

Difference_Q <- function(rmesh1, rmesh2, clean, normals, triangulate1, triangulate2) {
    .Call(`_Boov_Difference_Q`, rmesh1, rmesh2, clean, normals, triangulate1, triangulate2)
}

Union_EK <- function(rmeshes, clean, normals, triangulate) {
    .Call(`_Boov_Union_EK`, rmeshes, clean, normals, triangulate)
}

Union_Q <- function(rmeshes, clean, normals, triangulate) {
    .Call(`_Boov_Union_Q`, rmeshes, clean, normals, triangulate)
}

