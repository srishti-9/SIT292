
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