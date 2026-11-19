# Load required packages
library(glmnet)
library(boot)

# 1. Prepare predictor matrix (x) and response vector (y)
data(mtcars)
x <- as.matrix(mtcars[, -1])  # All columns except mpg
y <- mtcars$mpg

# 2. Find optimal lambda via cross-validation (alpha = 0 specifies Ridge)
set.seed(123)
cv_fit <- cv.glmnet(x, y, alpha = 0)
best_lambda <- cv_fit$lambda.min

# 3. Define the statistic function for boot()
ridge_boot_fn <- function(data, indices, lambda_val) {
  d <- data[indices, ] # Resample rows with replacement
  x_boot <- as.matrix(d[, -1])
  y_boot <- d$mpg
  
  # Fit ridge regression using the fixed optimal lambda
  fit <- glmnet(x_boot, y_boot, alpha = 0, lambda = lambda_val)
  
  # Return coefficients as a numeric vector
  return(as.vector(coef(fit)))
}

# 4. Run the bootstrap (e.g., R = 1000 resamples)
set.seed(123)
boot_results <- boot(
  data = mtcars,
  statistic = ridge_boot_fn,
  R = 1000,
  lambda_val = best_lambda
)

# 5. View bootstrap standard errors and 95% confidence intervals
print(boot_results)

# Example: 95% percentile confidence interval for the first predictor (cyl, index = 2)
boot.ci(boot_results, type = "perc", index = 2)
