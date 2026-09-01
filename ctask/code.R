# Module 9 - Using R for Linear Algebra
# Learning Evidence



# 1. Basic calculations: variables and arithmetic

a <- 3
b <- 7
the.value <- 12
the.index <- 2
the.index_2 <- 3

cat("a * b =", a * b, "\n")
cat("the.value * a =", the.value * a, "\n")
cat("a^the.index + b^the.index_2 =", a^the.index + b^the.index_2, "\n\n")

vec <- c(3, 2, 1, 0)
seq_range <- 4:8
seq_fn <- seq(4, 20, 5)

cat("vec:", vec, "\n")
cat("seq_range:", seq_range, "\n")
cat("seq_fn:", seq_fn, "\n\n")

x <- c(1, 6, 7, 9)
y <- c(-1, 2, 1, -2)

cat("x + y:", x + y, "\n")
cat("x - y:", x - y, "\n")
cat("x * y:", x * y, "\n")
cat("x / y:", x / y, "\n\n")



# 2. Data structures :  sample(), arrays, data frames, lists

my.vector <- c(1, 1, 1, 2, 3)
my.sample <- sample(my.vector, 20, replace = TRUE)

cat("Sample frequencies:\n")
print(table(my.sample))
cat("\n")

my.table <- array(0, c(3, 4))
my.table[1, ] <- 1:4
my.table[, 2] <- sample(c(-1, 1), 3, replace = TRUE)

cat("my.table after edits:\n")
print(my.table)
cat("\n")

my.data <- data.frame(
  x = 1:5,
  y = (1:5)^2,
  passed = sample(0:1, 5, replace = TRUE)
)

cat("Data frame:\n")
print(my.data)
cat("\n")

my.list <- list(
  a.vector = c(1, 3, 5, 7),
  a.matrix = matrix(c(1, 9, 5, 3, 2, 1), 2, 3),
  a.value = 42
)

cat("List element a.value:", my.list$a.value, "\n")
cat("List element a.matrix[2,2]:", my.list$a.matrix[2, 2], "\n\n")



# 3. Solving a 2x2 system of linear equations

# Two equations: ax + by = e, cx + dy = f
a <- 2; b <- 3; e <- 8      # equation 1
c <- 1; d <- -1; f <- 1     # equation 2

A <- matrix(c(a, c, b, d), nrow = 2)   # coefficient matrix
B <- c(e, f)                           # right-hand side

solution <- solve(A, B)
cat("Solution (x, y):\n")
print(solution)
cat("\n")

# Plot both lines to see where they cross
plot(0, 0, xlim = c(-10, 10), ylim = c(-10, 10), type = "n", xlab = "x", ylab = "y")
abline(a = e / b, b = -a / b, col = "red")
abline(a = f / d, b = -c / d, col = "blue")
points(solution[1], solution[2], pch = 19)



# 4. Gaussian elimination (3x3 system)

A <- matrix(c(2, 1, -1,
              -3, -1, 2,
              -2, 1, 2), nrow = 3, byrow = TRUE)
b <- c(8, -11, -3)
M <- cbind(A, b)

cat("Augmented matrix:\n")
print(M)
cat("\n")

M[2, ] <- M[2, ] - (M[2, 1] / M[1, 1]) * M[1, ]
M[3, ] <- M[3, ] - (M[3, 1] / M[1, 1]) * M[1, ]
M[3, ] <- M[3, ] - (M[3, 2] / M[2, 2]) * M[2, ]

cat("Triangular form:\n")
print(M)
cat("\n")

z <- M[3, 4] / M[3, 3]
y <- (M[2, 4] - M[2, 3] * z) / M[2, 2]
x <- (M[1, 4] - M[1, 2] * y - M[1, 3] * z) / M[1, 1]

cat("x =", x, " y =", y, " z =", z, "\n\n")



# 5. Matrix addition, scalar multiplication, matrix multiplication

A <- matrix(c(4, -2, 7, 1, 0, 5), nrow = 2, ncol = 3, byrow = TRUE)
B <- matrix(c(-3, 6, 2, 8, -1, 4), nrow = 2, ncol = 3, byrow = TRUE)

cat("A + B:\n")
print(A + B)
cat("\n")

k <- 3
cat("k * A:\n")
print(k * A)
cat("\n")

A <- matrix(c(0, 1, 0,
              1, 0, 0,
              0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)

B <- matrix(c(2, -1,
              4,  3,
              0,  5), nrow = 3, ncol = 2, byrow = TRUE)

cat("A %*% B:\n")
print(A %*% B)
cat("\n")

A <- matrix(c(1, 2,
              3, 4,
              5, 6), nrow = 3, ncol = 2, byrow = TRUE)

B <- matrix(c(1, 0, 2, 1,
              2, 1, 0, 3), nrow = 2, ncol = 4, byrow = TRUE)

cat("A %*% B (second example):\n")
print(A %*% B)
cat("\n")



# 6. Matrix operations summary via base R functions

A <- matrix(
  c(2, 1, 0,
    1, 3, 1,
    0, 1, 2),
  nrow = 3,
  byrow = TRUE
)

B <- matrix(c(1, 2, 3), nrow = 3)

cat("Transpose of A:\n")
print(t(A))
cat("\n")

cat("Determinant of A:", det(A), "\n\n")

cat("Inverse of A:\n")
print(solve(A))
cat("\n")

cat("Solving A X = B:\n")
print(solve(A, B))
cat("\n")

cat("Rank of A:", qr(A)$rank, "\n\n")



# 7. If/else, loops, Fibonacci, Collatz conjecture

n <- sample(1:10, 1)

if (n < 5) {
  cat(n, "is Bad\n")
} else if (n < 8) {
  cat(n, "is Okay\n")
} else {
  cat(n, "is Good\n")
}
cat("\n")

fib <- c(1, 1)

for (i in 1:10) {
  fib <- c(fib, fib[i] + fib[i + 1])
}

cat("Fibonacci sequence:\n")
print(fib)
cat("\n")

collatz <- function(start) {

  n <- start
  count <- 0

  while (n > 1) {

    if (n %% 2 == 0) {
      n <- n / 2
    } else {
      n <- 3 * n + 1
    }

    count <- count + 1
  }

  count
}

for (start_val in c(7, 27, 12, 19, 5)) {
  cat("Starting at", start_val, "took", collatz(start_val), "steps\n")
}
cat("\n")



# 8. Custom functions - handshakes, HCF

handshakes <- function(n) {
  n * (n - 1) / 2
}

cat("Handshakes among 6 people:", handshakes(6), "\n\n")

HCF <- function(a, b) {

  output <- 0
  div <- min(a, b)

  while (output == 0) {

    if (a %% div == 0) {
      if (b %% div == 0) {
        output <- div
      }
    }

    div <- div - 1
  }

  output
}

cat("HCF of 84 and 126:", HCF(84, 126), "\n\n")



# 9. Gram-Schmidt - custom function version (3x3, returns f1/f2/f3)
# No built-in matrix shortcuts used - only column references,
# sum() and sqrt()

gram_schmidt <- function(A) {

  a1 <- A[, 1]
  a2 <- A[, 2]
  a3 <- A[, 3]

  v1 <- a1

  v2 <- a2 -
    (sum(v1 * a2) / sum(v1 * v1)) * v1

  v3 <- a3 -
    (sum(v1 * a3) / sum(v1 * v1)) * v1 -
    (sum(v2 * a3) / sum(v2 * v2)) * v2

  f1 <- v1 / sqrt(sum(v1 * v1))
  f2 <- v2 / sqrt(sum(v2 * v2))
  f3 <- v3 / sqrt(sum(v3 * v3))

  list(f1 = f1, f2 = f2, f3 = f3)
}

A_test <- matrix(
  c(1, 1, 0,
    1, 0, 1,
    0, 1, 1),
  nrow = 3
)

gs_result <- gram_schmidt(A_test)

cat("f1:\n"); print(round(gs_result$f1, 4))
cat("f2:\n"); print(round(gs_result$f2, 4))
cat("f3:\n"); print(round(gs_result$f3, 4))
cat("\n")

Q_mine <- cbind(gs_result$f1, gs_result$f2, gs_result$f3)

cat("Checking Q^T Q (should be the identity):\n")
print(round(t(Q_mine) %*% Q_mine, 4))
cat("\n")



# 10. Gram-Schmidt - general n x n version (random matrix)

n <- 3

X <- matrix(sample(-5:5, n^2, replace = TRUE), nrow = n, ncol = n)

cat("Original vectors:\n")
print(X)
cat("\n")

if (abs(det(X)) < 1e-10) {
  stop(
    "The vectors are linearly dependent, so an orthonormal basis ",
    "cannot be generated from this set."
  )
}

Q <- matrix(0, nrow = n, ncol = n)

for (k in seq_len(n)) {

  current <- X[, k]

  if (k > 1) {
    for (j in seq_len(k - 1)) {
      projection <- sum(current * Q[, j])
      current <- current - projection * Q[, j]
    }
  }

  magnitude <- sqrt(sum(current^2))

  if (magnitude < 1e-10) {
    stop("A zero vector was produced during Gram-Schmidt.")
  }

  Q[, k] <- current / magnitude
}

cat("Final orthonormal matrix Q:\n")
print(round(Q, 4))
cat("\n")

cat("Checking Q^T Q:\n")
print(round(t(Q) %*% Q, 4))
cat("\n")



# 11. QR factorisation
# Reusing the custom Gram-Schmidt function: since A = QR and Q has
# orthonormal columns, t(Q) %*% A = t(Q) %*% Q %*% R = R

R <- t(Q_mine) %*% A_test

cat("R:\n")
print(round(R, 4))
cat("\n")

cat("Checking Q %*% R against the original matrix:\n")
print(round(Q_mine %*% R, 4))

if (isTRUE(all.equal(
  as.vector(Q_mine %*% R),
  as.vector(A_test),
  tolerance = 1e-8
))) {
  cat("Verified: Q %*% R = A\n")
}
cat("\n")



# 12. Gauss-Jordan elimination: inverse of A and solving Ax = b

n <- 3

repeat {
  M <- matrix(sample(-6:6, n^2, replace = TRUE), nrow = n)
  if (abs(det(M)) > 1e-10) break
}

x <- sample(-6:6, n, replace = TRUE)
b <- as.vector(M %*% x)

aug <- cbind(M, diag(n), b)

cat("Initial augmented matrix [A | I | b]\n\n")
print(aug)
cat("\n")

swap_rows <- function(mat, r1, r2) {
  temp <- mat[r1, ]
  mat[r1, ] <- mat[r2, ]
  mat[r2, ] <- temp
  mat
}

for (pivot in seq_len(n)) {

  if (abs(aug[pivot, pivot]) < 1e-10) {

    # guard against pivot being the last row, where (pivot+1):n
    # would otherwise count backwards instead of being empty
    if (pivot < n) {

      replacement <- which(
        abs(aug[(pivot + 1):n, pivot]) > 1e-10
      )

      if (length(replacement) == 0) {
        stop("Matrix is singular.")
      }

      new_row <- replacement[1] + pivot
      aug <- swap_rows(aug, pivot, new_row)

      cat("Swapped rows", pivot, "and", new_row, "\n")

    } else {
      stop("Matrix is singular.")
    }
  }

  if (pivot < n) {

    for (row in (pivot + 1):n) {

      if (abs(aug[row, pivot]) > 1e-10) {

        multiplier <- aug[row, pivot] / aug[pivot, pivot]
        aug[row, ] <- aug[row, ] - multiplier * aug[pivot, ]
      }
    }
  }
}

for (row in seq_len(n)) {

  pivot_value <- aug[row, row]

  if (abs(pivot_value) < 1e-10) {
    stop("Zero pivot encountered.")
  }

  aug[row, ] <- aug[row, ] / pivot_value
}

for (pivot in n:2) {
  for (row in (pivot - 1):1) {

    coefficient <- aug[row, pivot]

    if (abs(coefficient) > 1e-10) {
      aug[row, ] <- aug[row, ] - coefficient * aug[pivot, ]
    }
  }
}

inverse_M <- aug[, (n + 1):(2 * n)]
calculated_x <- aug[, 2 * n + 1]

cat("Inverse matrix:\n")
print(round(inverse_M, 3))

cat("\nSolution obtained from Gauss-Jordan elimination:\n")
print(round(calculated_x, 3))

cat("\nVerification of Ax = b:\n")
print(round(M %*% calculated_x, 3))
cat("\n")



# 13. Perpendicular vector to a given 2D vector

u <- sample(-6:6, 2, replace = TRUE)

cat("Given vector u = (", u[1], ", ", u[2], ")\n", sep = "")

v <- c(-u[2], u[1])

cat("A perpendicular vector is v = (",
    v[1], ", ", v[2], ")\n", sep = "")

dot_product <- sum(u * v)

cat("Dot product u . v =", dot_product, "\n")

if (dot_product == 0) {
  cat("Therefore, u and v are orthogonal.\n")
}
cat("\n")

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

abline(h = 0)
abline(v = 0)

arrows(0, 0, u[1], u[2], length = 0.1, col = "red", lwd = 2)
arrows(0, 0, v[1], v[2], length = 0.1, col = "blue", lwd = 2)



# 14. LU decomposition

size <- 3

A <- matrix(sample(-6:6, size^2, replace = TRUE), nrow = size)

while (abs(det(A)) < 1e-10) {
  A <- matrix(sample(-6:6, size^2, replace = TRUE), nrow = size)
}

L <- diag(size)
U <- matrix(0, nrow = size, ncol = size)

cat("Matrix A:\n")
print(A)
cat("\n")

for (k in 1:size) {

  for (j in k:size) {

    # seq_len(k - 1) is used instead of 1:(k - 1) so that when
    # k = 1 this correctly gives an empty sequence rather than
    # counting backwards as c(1, 0)
    U[k, j] <- A[k, j] -
      sum(L[k, seq_len(k - 1)] * U[seq_len(k - 1), j])
  }

  if (k < size) {

    for (i in (k + 1):size) {

      L[i, k] <- (
        A[i, k] -
        sum(L[i, seq_len(k - 1)] * U[seq_len(k - 1), k])
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
cat("\n")



# 15. Rank of a matrix via row reduction

A <- matrix(
  c(2, 3, 1, -1,
    3, 0, 2, 1,
    4, 6, 2, -2,
    1, -3, 1, 2),
  nrow = 4,
  byrow = TRUE
)

rows <- nrow(A)
cols <- ncol(A)
pivot_row <- 1

for (pivot_col in seq_len(cols)) {

  if (pivot_row > rows) break

  possible <- which(abs(A[pivot_row:rows, pivot_col]) > 1e-10)

  if (length(possible) == 0) next

  selected <- possible[1] + pivot_row - 1

  if (selected != pivot_row) {
    temporary <- A[pivot_row, ]
    A[pivot_row, ] <- A[selected, ]
    A[selected, ] <- temporary
  }

  if (pivot_row < rows) {
    for (r in (pivot_row + 1):rows) {

      multiplier <- A[r, pivot_col] / A[pivot_row, pivot_col]

      if (abs(multiplier) > 1e-10) {
        A[r, ] <- A[r, ] - multiplier * A[pivot_row, ]
      }
    }
  }

  pivot_row <- pivot_row + 1
}

nonzero_rows <- 0

for (r in seq_len(rows)) {
  if (any(abs(A[r, ]) > 1e-10)) {
    nonzero_rows <- nonzero_rows + 1
  }
}

cat("Row-echelon form:\n")
print(round(A, 3))
cat("\nRank =", nonzero_rows, "\n\n")



# 16. Power method - dominant eigenvalue/eigenvector

n <- 3

A_power <- matrix(sample(-5:5, n^2, replace = TRUE), nrow = n)

cat("Matrix for power method:\n")
print(A_power)
cat("\n")

x <- matrix(rep(1, n), nrow = n)

for (i in 1:25) {
  x <- A_power %*% x
  x <- x / sqrt(sum(x^2))
}

eigenvalue_estimate <- as.numeric(
  (t(x) %*% A_power %*% x) / (t(x) %*% x)
)

cat("Estimated dominant eigenvector:\n")
print(round(x, 4))
cat("\nEstimated dominant eigenvalue:", round(eigenvalue_estimate, 4), "\n\n")



# 17. QR algorithm - approximating all eigenvalues

A_qr <- A_power

for (i in 1:100) {

  decomposition <- qr(A_qr)
  Q_step <- qr.Q(decomposition)
  R_step <- qr.R(decomposition)

  A_qr <- R_step %*% Q_step
}

cat("Matrix after QR algorithm iterations (should be close to\n")
cat("upper triangular, with eigenvalues on the diagonal):\n")
print(round(A_qr, 4))

cat("\nDiagonal entries (eigenvalue estimates):\n")
print(round(diag(A_qr), 4))
cat("\n")



# 18. Eigenvalues and eigenvectors, with verification

n <- 3

A <- matrix(sample(-5:5, n^2, replace = TRUE), nrow = n, byrow = TRUE)

cat("Matrix A:\n")
print(A)
cat("\n")

result <- eigen(A)
values <- result$values
vectors <- result$vectors

cat("Eigenvalues:\n")
print(values)

cat("\nEigenvectors:\n")
print(vectors)
cat("\n")

cat("Verification:\n")

for (k in seq_len(n)) {

  lambda <- values[k]
  vector <- vectors[, k]

  left_side <- A %*% vector
  right_side <- lambda * vector

  cat("\nEigenvalue", k, "=", lambda, "\n")
  cat("A v =\n"); print(left_side)
  cat("lambda v =\n"); print(right_side)

  if (isTRUE(all.equal(
    as.vector(left_side),
    as.vector(right_side),
    tolerance = 1e-8
  ))) {
    cat("Verified: A v = lambda v\n")
  }
}