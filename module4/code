
# Determinant of a 3 x 3 matrix

size <- 3

mat <- matrix(
  sample(-9:9, size^2, replace = TRUE),
  nrow = size
)

cat("Matrix:\n")
print(mat)
cat("\n")

# Calculate determinant directly
answer <- det(mat)

cat("Determinant =", answer, "\n")


#LU decomposition

size <- 3

A <- matrix(
  sample(-6:6, size^2, replace = TRUE),
  nrow = size
)

while (abs(det(A)) < 1e-10) {
  A <- matrix(
    sample(-6:6, size^2, replace = TRUE),
    nrow = size
  )
}

L <- diag(size)
U <- matrix(0, nrow = size, ncol = size)

cat("Matrix A:\n")
print(A)
cat("\n")

for (k in 1:size) {

  # Calculate the kth row of U
  for (j in k:size) {

    U[k, j] <- A[k, j] -
      sum(L[k, 1:(k - 1)] * U[1:(k - 1), j])
  }

  # Calculate the kth column of L
  if (k < size) {

    for (i in (k + 1):size) {

      L[i, k] <- (
        A[i, k] -
        sum(L[i, 1:(k - 1)] * U[1:(k - 1), k])
      ) / U[k, k]
    }
  }
}

cat("L matrix:\n")
print(round(L, 4))

cat("\nU matrix:\n")
print(round(U, 4))

cat("\nChecking whether A = LU:\n")
print(round(L %*% U, 4))