
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

