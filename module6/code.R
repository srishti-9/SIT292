
# Finding the rank of a matrix using row reduction

A <- matrix(
  c(2, 3, 1, -1,
    3, 0, 2, 1,
    4, 6, 2, -2,
    1, -3, 1, 2),
  nrow = 4,
  byrow = TRUE
)

cat("Original matrix:\n")
print(A)
cat("\n")

rows <- nrow(A)
cols <- ncol(A)
pivot_row <- 1

# Gaussian elimination
for (pivot_col in seq_len(cols)) {

  if (pivot_row > rows) {
    break
  }

  # Find a non-zero entry that can be used as the pivot
  possible <- which(
    abs(A[pivot_row:rows, pivot_col]) > 1e-10
  )

  if (length(possible) == 0) {
    next
  }

  selected <- possible[1] + pivot_row - 1

  # Exchange rows when required
  if (selected != pivot_row) {
    temporary <- A[pivot_row, ]
    A[pivot_row, ] <- A[selected, ]
    A[selected, ] <- temporary

    cat("Rows", pivot_row, "and", selected, "exchanged.\n")
  }

  # Eliminate entries below the pivot
  if (pivot_row < rows) {
    for (r in (pivot_row + 1):rows) {

      multiplier <- A[r, pivot_col] / A[pivot_row, pivot_col]

      if (abs(multiplier) > 1e-10) {
        A[r, ] <- A[r, ] -
          multiplier * A[pivot_row, ]
      }
    }
  }

  pivot_row <- pivot_row + 1
}

cat("Row-echelon form:\n")
print(round(A, 3))
cat("\n")

# Count non-zero rows
nonzero_rows <- 0

for (r in seq_len(rows)) {
  if (any(abs(A[r, ]) > 1e-10)) {
    nonzero_rows <- nonzero_rows + 1
  }
}

cat("Rank =", nonzero_rows, "\n")

# Eigenvalues and eigenvectors

n <- 3

A <- matrix(
  sample(-5:5, n^2, replace = TRUE),
  nrow = n,
  byrow = TRUE
)

cat("Matrix A:\n")
print(A)
cat("\n")

# Calculate eigenvalues and eigenvectors
result <- eigen(A)

values <- result$values
vectors <- result$vectors

cat("Eigenvalues:\n")
print(values)

cat("\nEigenvectors:\n")
print(vectors)


#Verify Av = lambda*v


cat("\nVerification:\n")

for (k in seq_len(n)) {

  lambda <- values[k]
  vector <- vectors[, k]

  left_side <- A %*% vector
  right_side <- lambda * vector

  cat("\nEigenvalue", k, "=", lambda, "\n")

  cat("A v =\n")
  print(left_side)

  cat("lambda v =\n")
  print(right_side)

  if (isTRUE(all.equal(
    as.vector(left_side),
    as.vector(right_side),
    tolerance = 1e-8
  ))) {
    cat("Verified: A v = lambda v\n")
  }
}