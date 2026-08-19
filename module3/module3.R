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