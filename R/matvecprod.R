#' Title
#'
#' @param A 
#' @param B 
#' @param V 
#' @param W 
#'
#' @return
#' @export
#'
#' @examples
matvecprod <- function (A, B, V, W) {
  m = dim(A)[1]
  p = dim(B)[2]
  r = dim(B)[1]
  if (W == 1){
  if (is.matrix(A) & is.matrix(B) == F){
    warning("Error! The first two arguments must be n by m matrices")
  } else {
    if (dim(A)[2] != dim(B)[1]){
    warning("Error! Check that the number of columns and rows match in your matrix")} else {
        C_1 = matrix(data = rep(0, m*p), nrow=m, ncol=p)
      for (i in 1:m)
      for (j in 1:p){
            C_1[i, j] = sum(A[i, ] * B[, j])
          }
      if (is.vector(V) == F){
    warning("Error! Third Argument must be a Vector")
      } else {
        n = length(V)
        if(n != dim(C_1)[2]){
          warning("Error! Check number of columns and dimensions of vector")} else {
              M_1 = rep(0, m)
            for (i in 1:m){
                M_1[i] = sum(C_1[i, ] * V)
              }
            }
          }
        }
  }
  return(M_1)
  }
  if (W == 2){
    if (is.matrix(A) & is.matrix(B) == F){
      warning("Error! The first two arguments must be n by m matrices")
    } else {
    if (is.vector(V) == F){
      warning("Error! The third argument must be a vector")
    } else {
      n = length(V)
      if(n != p){
        warning("Error! Check number of columns and dimensions of vector")} else {
          C_2 = rep(0, r)
          for (i in 1:r){
            C_2[i] = sum(B[i, ] * V)
          }
          if (dim(A)[2] != length(C_2)){
            warning("Error! Check number of columns and dimensions of vector")} else {
              M_2 = rep(0, m)
              for (i in 1:m){
                M_2[i] = sum(A[i, ] * C_2)
              }
            }
          return(M_2)
    }
    }
  }
  }
}


    



