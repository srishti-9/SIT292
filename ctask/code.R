
## SIT292  Using R for Linear Algebra
## Learning Evidence

## 1. Solving a 2x2 system of linear equations 

# Two equations: ax + by = e, cx + dy = f
a <- 2; b <- 3; e <- 8      # equation 1
c <- 1; d <- -1; f <- 1     # equation 2

# Solve using matrices: A %*% x = B  ->  x = solve(A) %*% B
A <- matrix(c(a, c, b, d), nrow = 2)   # coefficient matrix
B <- c(e, f)                           # right-hand side

solution <- solve(A, B)   # solves for x and y
print(solution)           # solution[1] = x, solution[2] = y

# Plot both lines to see where they cross
plot(0, 0, xlim = c(-10,10), ylim = c(-10,10), type = "n", xlab = "x", ylab = "y")
abline(a = e/b, b = -a/b, col = "red")       # line 1
abline(a = f/d, b = -c/d, col = "blue")      # line 2
points(solution[1], solution[2], pch = 19)   # mark the solution point


## 2. Gaussian elimination (3x3 system) 

#gaussian elimination
# Augmented matrix [A | b]
A <- matrix(c(2, 1, -1,
              -3, -1, 2,
              -2, 1, 2), nrow = 3, byrow = TRUE)
b <- c(8, -11, -3)
M <- cbind(A, b)   # combine into one matrix
print(M)

# Step 1: clear column 1 below row 1
M[2,] <- M[2,] - (M[2,1]/M[1,1]) * M[1,]
M[3,] <- M[3,] - (M[3,1]/M[1,1]) * M[1,]

# Step 2: clear column 2 below row 2
M[3,] <- M[3,] - (M[3,2]/M[2,2]) * M[2,]

print(M)   # now in triangular form

# Step 3: back-substitution to solve
z <- M[3,4] / M[3,3]
y <- (M[2,4] - M[2,3]*z) / M[2,2]
x <- (M[1,4] - M[1,2]*y - M[1,3]*z) / M[1,1]

cat("x =", x, " y =", y, " z =", z, "\n")


##  3. Matrix addition 

# Matrix Addition (A + B)
A <- matrix(c(4, -2, 7, 1, 0, 5), nrow = 2, ncol = 3, byrow = TRUE)
B <- matrix(c(-3, 6, 2, 8, -1, 4), nrow = 2, ncol = 3, byrow = TRUE)

print("A =")
print(A)
print("B =")
print(B)
print("A + B =")
print(A + B)


##  4. Scalar multiplication 

# Scalar Multiplication (k * A)
k <- 3
A <- matrix(c(1, -5, 2, 4, 0, -3), nrow = 2, ncol = 3, byrow = TRUE)

print("A =")
print(A)
print(paste("k =", k))
print("k * A =")
print(k * A)


##  5. Matrix multiplication 

#  Matrix Multiplication (A x B)
A <- matrix(c(0, 1, 0,
              1, 0, 0,
              0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)

B <- matrix(c(2, -1,
              4,  3,
              0,  5), nrow = 3, ncol = 2, byrow = TRUE)

print("A =")
print(A)
print("B =")
print(B)
print("A x B =")
print(A %*% B)


##  6. Matrix multiplication (example 2) 

A <- matrix(c(1, 2,
              3, 4,
              5, 6), nrow = 3, ncol = 2, byrow = TRUE)

B <- matrix(c(1, 0, 2, 1,
              2, 1, 0, 3), nrow = 2, ncol = 4, byrow = TRUE)

print("A =")
print(A)
print("B =")
print(B)
print("A x B =")
print(A %*% B)


##  7. Gauss-Jordan elimination: inverse + solving Ax = b 

# Matrix size
n <- 3

repeat {
  M <- matrix(
    sample(-6:6, n^2, replace = TRUE),
    nrow = n
  )

  if (abs(det(M)) > 1e-10) {
    break
  }
}

# Create a random solution vector
x <- sample(-6:6, n, replace = TRUE)

# Calculate b from Ax = b
b <- as.vector(M %*% x)

# Form the augmented matrix [A | I | b]
aug <- cbind(M, diag(n), b)

cat("Initial augmented matrix [A | I | b]\n\n")
print(aug)
cat("\n")


# Helper function for swapping two row

swap_rows <- function(mat, r1, r2) {
  temp <- mat[r1, ]
  mat[r1, ] <- mat[r2, ]
  mat[r2, ] <- temp
  return(mat)
}


# Gaussian elimination

for (pivot in seq_len(n)) {

  # If the pivot is zero, look for a suitable row below it
  if (abs(aug[pivot, pivot]) < 1e-10) {

    replacement <- which(
      abs(aug[(pivot + 1):n, pivot]) > 1e-10
    )

    if (length(replacement) == 0) {
      stop("Matrix is singular.")
    }

    new_row <- replacement[1] + pivot
    aug <- swap_rows(aug, pivot, new_row)

    cat("Swapped rows", pivot, "and", new_row, "\n")
    print(round(aug, 3))
    cat("\n")
  }

  if (pivot < n) {

    for (row in (pivot + 1):n) {

      if (abs(aug[row, pivot]) > 1e-10) {

        multiplier <- aug[row, pivot] / aug[pivot, pivot]

        aug[row, ] <- aug[row, ] -
          multiplier * aug[pivot, ]

        cat(
          "R", row, " <- R", row,
          " - (", round(multiplier, 3), ")R", pivot,
          "\n",
          sep = ""
        )
      }
    }

    print(round(aug, 3))
    cat("\n")
  }
}


# Convert the diagonal entries into 1

for (row in seq_len(n)) {

  pivot_value <- aug[row, row]

  if (abs(pivot_value) < 1e-10) {
    stop("Zero pivot encountered.")
  }

  aug[row, ] <- aug[row, ] / pivot_value
}

cat("After normalising the pivots:\n")
print(round(aug, 3))
cat("\n")


# Gauss-Jordan elimination
# Remove entries above each pivot

for (pivot in n:2) {

  for (row in (pivot - 1):1) {

    coefficient <- aug[row, pivot]

    if (abs(coefficient) > 1e-10) {

      aug[row, ] <- aug[row, ] -
        coefficient * aug[pivot, ]

      cat(
        "R", row, " <- R", row,
        " - (", round(coefficient, 3), ")R", pivot,
        "\n",
        sep = ""
      )
    }
  }

  print(round(aug, 3))
  cat("\n")
}

# Extract the results

inverse_M <- aug[, (n + 1):(2 * n)]
calculated_x <- aug[, 2 * n + 1]

cat("Inverse matrix:\n")
print(round(inverse_M, 3))

cat("\nOriginal solution vector:\n")
print(x)

cat("\nSolution obtained from Gauss-Jordan elimination:\n")
print(round(calculated_x, 3))

cat("\nVerification of Ax = b:\n")
print(round(M %*% calculated_x, 3))

cat("\nOriginal b:\n")
print(b)


##  8. Perpendicular vector to a given 2D vector 

# Finding a vector perpendicular to a given 2D vector

u <- sample(-6:6, 2, replace = TRUE)

cat("Given vector u = (", u[1], ", ", u[2], ")\n", sep = "")

# Rotate the vector by 90 degrees
v <- c(-u[2], u[1])

cat("A perpendicular vector is v = (",
    v[1], ", ", v[2], ")\n", sep = "")

# Check that the dot product is zero
dot_product <- sum(u * v)

cat("Dot product u . v =", dot_product, "\n")

if (dot_product == 0) {
  cat("Therefore, u and v are orthogonal.\n")
}



# Plot the two vectors

# Determine suitable plot limits
maximum <- max(abs(c(u, v))) + 2

plot(
  0, 0,
  type = "n",
  xlim = c(-maximum, maximum),
  ylim = c(-maximum, maximum),
  xlab = "x",
  ylab = "y",
  asp = 1,
  main = "Orthogonal Vectors"
)

# Coordinate axes
abline(h = 0)
abline(v = 0)

# Original vector
arrows(
  0, 0,
  u[1], u[2],
  length = 0.1,
  col = "red",
  lwd = 2
)

# Perpendicular vector
arrows(
  0, 0,
  v[1], v[2],
  length = 0.1,
  col = "blue",
  lwd = 2
)


##  9. Determinant of a 3x3 matrix 

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


##  10. LU decomposition 

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


##  11. Gram-Schmidt orthonormalisation 

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


##  12. Rank of a matrix via row reduction 

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


##  13. Eigenvalues and eigenvectors 

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