library(tidyverse)

diabetes <- read.csv("http://hastie.su.domains/CASI_files/DATA/diabetes.csv")

#' The variable names are:
#' 
names(diabetes)
str(diabetes)

# Diabetes data of Section 7.3
# These data consist of observations on 442 patients, with the response of interest being 
# a quantitative measure of disease progression one year after baseline: prog
# There are ten baseline variables---age, sex, body-mass index, average blood pressure, and six blood serum measurements.
# 
# First used in LARS paper
# 
# Please note: in table 7.2, we standardized the centered predictor variables to be unit L2 norm. 
# In table 20.1 we standardized to be unit variance. The unit norm standardization was inherited from 
# our work on LARS, where Euclidean geometry played a role. 

diab <- diabetes |>
  select(-X) |>
  relocate(prog, .before = age)
str(diab)

save(diab, file = "diab.RData")

use_data_doc(diab, file="diab.Rd")

load("data/diab.RData")

# standardize X to mean zero, sum of squares = 1
colMeans(diab[, -1])

diab_scaled <- scale(diab) |> as.data.frame()
colSums(diab_scaled^2)

diab_scaled <- diab
diab_scaled[, -1] <- as.data.frame(lapply(diab[, -1], function(x) x - mean(x)))

scale2 <- function(x, na.rm = FALSE){ 
  x <- (x - mean(x, na.rm = na.rm))
  x <- x / sum(x^2)
}

diab_scaled <- diab |>
  mutate(across(!"prog"), scale2())

diab.lm <- lm(prog ~ ., data = diab_scaled)
coef(diab.lm)


lambda <- c(0, 0.005, 0.01, 0.02, 0.05, 0.075, 0.1, 0.15)
diab.lridge <- ridge(prog ~ ., data = diab, lambda =lambda)

