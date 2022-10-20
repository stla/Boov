isPositiveInteger <- function(x){
  is.numeric(x) && length(x) == 1L && !is.na(x) && floor(x) == x
}

isStrictPositiveInteger <- function(x){
  isPositiveInteger(x) && x > 0
}

#' @importFrom gmp `%*%`
#' @noRd
`%^%` <- function(A, n){
  Reduce(gmp::`%*%`, replicate(n, A, simplify = FALSE))
}

#' @title Rational approximation of square roots
#' @description Returns a rational approximation of the square root of
#'   an integer.
#'
#' @param x the positive integer whose square root is desired
#' @param n a positive integer, the higher the better approximation
#' @param ... ignored
#'
#' @return The \code{qsqrt} function returns a \strong{gmp} rational number
#'   (class \code{\link[gmp]{bigq}}) approximating the square root of
#'   \code{x}. The \code{qsqrt2}, \code{qsqrt3}, and \code{qsqrtPhi} functions
#'   return a \strong{gmp} rational number approximating the square root of
#'   \code{2}, \code{3}, and \code{phi} (the golden number) respectively.
#'   Their value converge more fastly than the value obtained with \code{qsqrt}.
#'
#' @importFrom gmp as.bigz as.bigq asNumeric matrix.bigq add.bigq factorialZ `%*%`
#'
#' @aliases qsqrt qsqrt2 qsqrt3 qsqrtPhi print.qsqrt
#' @rdname qsqrt
#' @export
#'
#' @examples
#' library(Boov)
#' qsqrt(2, 7)
#' qsqrt2(7)
#' qsqrt3(22)
#' qsqrtPhi(17)
qsqrt <- function(x, n){
  stopifnot(isPositiveInteger(x))
  stopifnot(isStrictPositiveInteger(n))
  zero <- as.bigz(0L)
  one <- as.bigz(1L)
  A <- matrix.bigq(
    c(zero, as.bigz(x)-1L, one, as.bigz(2L)),
    nrow = 2L, ncol = 2L
  )
  zs <- c(gmp::`%*%`(A %^% n, c(zero, one)))
  out <- as.bigq(zs[2L], zs[1L]) - 1L
  attr(out, "error") <- abs(asNumeric(out) - sqrt(x))
  class(out) <- c("qsqrt", class(out))
  out
}

#' @rdname qsqrt
#' @export
qsqrt2 <- function(n){
  stopifnot(isStrictPositiveInteger(n) && n >= 2)
  out <- as.bigq(99L, 70L) + Reduce(add.bigq, sapply(2L:n, function(i){
    numer <- 10L * prod(1L - 2L*(0L:(i-1L)))
    denom <- 7L * factorialZ(i) * as.bigz(-100L)^i
    as.bigq(numer, denom)
  }))
  attr(out, "error") <- abs(asNumeric(out) - sqrt(2))
  class(out) <- c("qsqrt", class(out))
  out
}

#' @rdname qsqrt
#' @export
qsqrt3 <- function(n){
  stopifnot(isStrictPositiveInteger(n) && n >= 2)
  out <- as.bigq(7L, 4L) + Reduce(add.bigq, sapply(2L:n, function(i){
    as.bigq(
      2L * prod(1L - 2L*(0L:(i-1L))),
      factorialZ(i) * as.bigz(-8L)^i
    )
  }))
  attr(out, "error") <- abs(asNumeric(out) - sqrt(3))
  class(out) <- c("qsqrt", class(out))
  out
}

#' @rdname qsqrt
#' @export
qsqrtPhi <- function(n){
  stopifnot(isStrictPositiveInteger(n) && n >= 2)
  out <- as.bigq(13L, 8L) + Reduce(add.bigq, sapply(2L:n, function(i){
    as.bigq(
      10L * prod(1L - 2L*(0L:(i-1L))),
      as.bigz(8L) * factorialZ(i) * as.bigz(-10L)^(i)
    )
  }))
  attr(out, "error") <- abs(asNumeric(out) - (1+sqrt(5))/2)
  class(out) <- c("qsqrt", class(out))
  out
}

#' @rdname qsqrt
#' @exportS3Method print qsqrt
print.qsqrt <- function(x, ...){
  print(as.bigq(x))
  cat('attr("error")')
  print(attr(x, "error"))
  invisible(NULL)
}
