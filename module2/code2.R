# Matrix Addition (A + B) 
A <- matrix(c(4, -2, 7, 1, 0, 5), nrow = 2, ncol = 3, byrow = TRUE)
B <- matrix(c(-3, 6, 2, 8, -1, 4), nrow = 2, ncol = 3, byrow = TRUE)

print("A =")
print(A)
print("B =")
print(B)
print("A + B =")
print(A + B)


# Scalar Multiplication (k * A) 
k <- 3
A <- matrix(c(1, -5, 2, 4, 0, -3), nrow = 2, ncol = 3, byrow = TRUE)

print("A =")
print(A)
print(paste("k =", k))
print("k * A =")
print(k * A)


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


#  Another Matrix Multiplication Example ---
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