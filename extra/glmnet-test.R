library(glmnet)

longley.y <- longley[, "Employed"]
longley.X <- data.matrix(longley[, c(2:6,1)])

lambda <- c(0, 0.001, 0.005, 0.01, 0.02, 0.04, 0.08)
lambdaf <- c(expression(~widehat(beta)^OLS), lambda[-1])
gridge <- glmnet(longley.X, longley.y, alpha=0, lambda=lambda)
gridge

coef(gridge)

plot(gridge)
