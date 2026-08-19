
# Gram-Schmidt orthonormalisation

n <- 3

X <- matrix(
  sample(-5:5, n^2, replace = TRUE),
  nrow = n,
  ncol = n
)

cat("Original vectors:\n")
print(X)
cat("\n")



# Check that the vectors are linearly independent

if (abs(det(X)) < 1e-10) {

  stop(
    "The vectors are linearly dependent, so an orthonormal basis ",
    "cannot be generated from this set."
  )
}



# Gram-Schmidt process


Q <- matrix(0, nrow = n, ncol = n)

for (k in seq_len(n)) {

  # Begin with the current vector
  current <- X[, k]

  # Remove the components in the directions
  # of all previously calculated orthonormal vectors
  if (k > 1) {

    for (j in seq_len(k - 1)) {

      projection <- sum(current * Q[, j])

      current <- current -
        projection * Q[, j]
    }
  }

  # Calculate its magnitude
  magnitude <- sqrt(sum(current^2))

  if (magnitude < 1e-10) {
    stop("A zero vector was produced during Gram-Schmidt.")
  }

  # Normalise the vector
  Q[, k] <- current / magnitude

  cat("Orthonormal vector", k, ":\n")
  print(round(Q[, k], 4))
  cat("\n")
}



# Display the resulting orthonormal set

cat("Final orthonormal matrix Q:\n")
print(round(Q, 4))
cat("\n")



# Verification


cat("Checking Q^T Q:\n")
print(round(t(Q) %*% Q, 4))

cat("\nThe diagonal entries should be 1 and\n")
cat("the off-diagonal entries should be 0.\n")