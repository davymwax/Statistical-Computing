matvecprod <- function (A, B, V) {
  if (is.matrix(A) & is.matrix(B) == F){
    warning("Error! Arguments must be n by m matrices")
  } else {
    if (dim(A)[2] != dim(B)[1]){
    warning("Error! Matrix Operation Undefined")} else {
        m = dim(A)[1]
        p = dim(B)[2]
        C = matrix(data = rep(0, m*p), nrow=m, ncol=p)
      for (i in 1:m)
      for (j in 1:p){
            C[i, j] = sum(A[i, ] * B[, j])
          }
      if (is.vector(V) == F){
    warning("Error! Third Argument must be a Vector")
      } else {
        p = dim(C)[2]
        n = length(V)
        if(n != p){
          warning("Error! Operation Undefined")} else {
              M = rep(0, m)
            for (i in 1:m){
                M[i] = sum(C[i, ] * V)
              }
            }
          }
        }
  }
  return(M)
  }
  

matvecprod(A, B, V)
    



