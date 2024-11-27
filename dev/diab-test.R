# -------------------------------------

library(tidyverse)
library(genridge)
library(lmridge)


load(here::here("data", "diab.RData"))

# standardize X to mean zero, sum of squares = 1


scale2 <- function(x, na.rm = FALSE){ 
  x <- (x - mean(x, na.rm = na.rm))
  x <- x / sqrt(sum(x^2))
  x
}

xt <- 1:10
st <- scale2(xt)
sum(st^2)

diab_scaled <- diab |>
  mutate(across(age:glu, scale2))

diab.lm <- lm(prog ~ ., data = diab_scaled)
coef(diab.lm)


lambda <- c(0, 0.005, 0.01, 0.02, 0.05, 0.075, 0.1, 0.15)
diab.ridge <- ridge(prog ~ ., data = diab_scaled, lambda =lambda)
diab.ridge

# this gives CASI results, Table 7.3
diab.lmridge <- lmridge(prog ~ ., data = diab_scaled, K =lambda, scaling="sc")
plot(diab.lmridge)

#lmridge(prog ~ ., data = diab_scaled, K =lambda, scaling="scaled")
