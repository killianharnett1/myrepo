library(IMIFA)

# Q3.1(a)
data("olive", package = "IMIFA")
original_class <- class(olive)
cat("Original class of olive dataset:", original_class, "\n")
new_class <- c("mydf", original_class)
class(olive) <- new_class
cat("Updated class of olive dataset:", class(olive), "\n")

# Q3.1(b)
print.mydf <- function(x, n = 5, ...) {
  if (!inherits(x, "mydf")) stop("Object is not of class 'mydf'")
  rows <- min(n, nrow(x))
  cols <- min(n, ncol(x))
  if (n > nrow(x)) warning("The requested number of rows (n) exceeds the number of rows in the dataset.")
  if (n > ncol(x)) warning("The requested number of columns (n) exceeds the number of columns in the dataset.")
  cat(sprintf("Printing the first %d rows and %d columns of the data frame:\n", rows, cols))
  sub_data <- x[1:rows, 1:cols]
  print(as.data.frame(sub_data), ...)
  if (nrow(x) > rows) cat("... and", nrow(x) - rows, "more rows\n")
  if (ncol(x) > cols) cat("... and", ncol(x) - cols, "more columns\n")
}
print(olive)

# Q3.1(c)
print(olive[-1, ])

# Q3.1(d)
print.mydf <- function(x, n = 5, ...) {
  if (!inherits(x, "mydf")) stop("Object is not of class 'mydf'")
  rows <- min(n, nrow(x))
  cols <- min(n, ncol(x))
  if (n > nrow(x)) warning("The requested number of rows (n) exceeds the number of rows in the dataset.")
  if (n > ncol(x)) warning("The requested number of columns (n) exceeds the number of columns in the dataset.")
  cat(sprintf("Printing the first %d rows and %d columns of the data frame:\n", rows, cols))
  sub_data <- x[1:rows, 1:cols]
  print(as.data.frame(sub_data), ...)
  if (nrow(x) > rows) cat("... and", nrow(x) - rows, "more rows\n")
  if (ncol(x) > cols) cat("... and", ncol(x) - cols, "more columns\n")
}
print(olive, n = 1000)

# Q3.1(e)
max_dens <- function(obj, ...) UseMethod("max_dens")
max_dens.density <- function(obj, ...) {
  if (!inherits(obj, "density")) stop("The object must be of class 'density'")
  max_index <- which.max(obj$y)
  c(y = obj$y[max_index], x = obj$x[max_index])
}
set.seed(123)
data <- rnorm(100)
density_obj <- density(data)
max_peak <- max_dens(density_obj)
print(max_peak)

# Q3.1(f)
numeric_vars <- olive[, sapply(olive, is.numeric)]
density_peaks <- vapply(
  numeric_vars,
  function(var) {
    density_obj <- density(var, kernel = "epanechnikov")
    max_dens(density_obj)
  },
  numeric(2)
)
rownames(density_peaks) <- c("y (height)", "x (location)")
density_peaks

# Q3.1(g)
palmitic_density <- density(olive$palmitic, kernel = "epanechnikov")
max_peak <- max_dens(palmitic_density)
plot(palmitic_density, main = "Density Plot of Palmitic", xlab = "Palmitic", ylab = "Density", col = "black", lwd = 1)
segments(x0 = max_peak["x"], y0 = 0, x1 = max_peak["x"], y1 = max_peak["y"], col = "red", lwd = 1)

# Q3.2(a)
f <- function(x) x + 3
g <- function(x) x^2 - 4

# Q3.2(b)
curve(f, from = -4, to = 4, ylim = c(-6, 12), col = "blue", lwd = 2, ylab = "y", xlab = "x", main = "Visualization of f(x) and g(x)")
curve(g, from = -4, to = 4, col = "red", lwd = 2, add = TRUE)
legend("topright", legend = c("f(x) = x + 3", "g(x) = x^2 - 4"), col = c("blue", "red"), lwd = 2)

# Q3.2(c)
intersection_function <- function(x) f(x) - g(x)
root1 <- uniroot(intersection_function, interval = c(-4, 0), tol = 1e-08)$root
root2 <- uniroot(intersection_function, interval = c(0, 4), tol = 1e-08)$root
root1
root2

# Q3.2(d)
area1 <- integrate(function(x) g(x) - f(x), lower = root1, upper = 0)$value
area2 <- integrate(function(x) f(x) - g(x), lower = 0, upper = root2)$value
total_area <- area1 + area2
total_area

# Q3.3(a)
df <- with(mtcars, data.frame(y = mpg, x1 = disp, x2 = hp, x3 = wt))
head(df)

# Q3.3(b)
nll_lm <- function(data, par) {
  y <- data$y
  X <- as.matrix(cbind(1, data[, -1]))
  beta <- par[-length(par)]
  sigma <- par[length(par)]
  residuals <- y - X %*% beta
  if (sigma <= 0) return(Inf)
  -sum(dnorm(residuals, mean = 0, sd = sigma, log = TRUE))
}
initial_par <- c(0, 0, 0, 0, 1)
nll_value <- nll_lm(data = df, par = initial_par)
nll_value

# Q3.3(c)
optim_result <- optim(
  par = initial_par,
  fn = nll_lm,
  data = df,
  method = "L-BFGS-B",
  lower = c(rep(-Inf, ncol(df)), 0.001),
  upper = c(rep(Inf, ncol(df)), Inf)
)
optim_result$par

# Q3.3(e)
y <- df$y
X <- as.matrix(cbind(1, df[, -1]))
beta_matrix <- solve(t(X) %*% X) %*% t(X) %*% y
beta_optim <- optim_result$par[1:ncol(X)]
comparison <- data.frame(
  Method = c("Matrix", "Optim"),
  Intercept = c(beta_matrix[1], beta_optim[1]),
  x1 = c(beta_matrix[2], beta_optim[2]),
  x2 = c(beta_matrix[3], beta_optim[3]),
  x3 = c(beta_matrix[4], beta_optim[4])
)
comparison

# Q3.3(f)
residuals_matrix <- y - X %*% beta_matrix
RSS <- sum(residuals_matrix^2)
n <- nrow(X)
p <- ncol(X)
sigma_matrix <- sqrt(RSS / (n - p))
sigma_optim <- optim_result$par[length(optim_result$par)]
comparison_sigma <- data.frame(
  Method = c("Matrix", "Optim"),
  Sigma = c(sigma_matrix, sigma_optim)
)
comparison_sigma

# Q3.3(h)
optim_result <- optim(
  par = initial_par,
  fn = nll_lm,
  data = df,
  method = "L-BFGS-B",
  lower = c(rep(-Inf, ncol(df)), 0.001),
  upper = c(rep(Inf, ncol(df)), Inf),
  hessian = TRUE
)
cov_matrix <- solve(optim_result$hessian)
std_errors <- sqrt(diag(cov_matrix))
names(std_errors) <- c("Intercept", colnames(df)[2:ncol(df)], "Sigma")
std_errors

#Q3.4
lm_fit <- lm(y ~ x1 + x2 +x3, data=df)
lm_coefficients <- coef(lm_fit)
lm_sigma <- summary(lm_fit)$sigma

lm_coefficients
lm_sigma