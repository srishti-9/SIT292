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

