#' ---
#' title: Longley data example for the genridge package
#' ---

library(genridge)
library(car)     # for vif() generic

#' ## Load the data
data(longley, package="datasets")
str(longley)

#' ## Fit the model predicting `Employed` & examine VIFs
longley.lm <- lm(Employed ~ GNP + Unemployed + Armed.Forces + 
                   Population + Year + GNP.deflator, 
                 data=longley)
vif(longley.lm)

#' ## Get ridge regression estimates

lambda <- c(0, 0.002, 0.005, 0.01, 0.02, 0.04, 0.08)
lridge <- ridge(Employed ~ GNP + Unemployed + Armed.Forces + 
                  Population + Year + GNP.deflator, 
                data=longley, lambda=lambda)
print(lridge, digits = 2)

#' ## The `"ridge"` object

names(lridge)

#' ## Standard, univariate trace plot
traceplot(lridge, 
          X = "lambda",
          xlab = "Ridge constant (k)",
          xlim = c(-0.02, 0.08), 
          cex.lab=1.25)

#' ## Better: plot against effective degrees of freedom
traceplot(lridge, 
          X = "df",
          xlim = c(4, 6.2), cex.lab=1.25)

#' ## Even BETTER: Bivariate ridge trace plot
#' Shows bias _and_ precision by confidence ellipses for the coefficients

clr <-  c("black", "red", "brown", "darkgreen","blue", "cyan4", "magenta")
pch <- c(15:18, 7, 9, 12)
lambdaf <- c(expression(~widehat(beta)^OLS), as.character(lambda[-1]))

#' GNP vs. Unemployed
plot(lridge, variables=c(1,2),
     radius=0.5, cex.lab=1.5, col=clr, 
     labels=NULL, fill=TRUE, fill.alpha=0.2)
text(lridge$coef[1,1], lridge$coef[1,2], 
     expression(~widehat(beta)^OLS), cex=1.5, pos=4, offset=.1)
text(lridge$coef[-1,c(1,2)], lambdaf[-1], pos=3, cex=1.3)

#' GNP vs. Population
plot(lridge, variables=c(1,4),
     radius=0.5, cex.lab=1.5, col=clr, 
     labels=NULL, fill=TRUE, fill.alpha=0.2)
text(lridge$coef[1,1], lridge$coef[1,4], 
     expression(~widehat(beta)^OLS), cex=1.5, pos=4, offset=.1)
text(lridge$coef[-1,c(1,4)], lambdaf[-1], pos=3, cex=1.3)

#' ## Scatterplot matrix
pairs(lridge, radius=0.5, diag.cex = 2, 
      fill = TRUE, fill.alpha = 0.1)

#' ## Visualizing the bias-variance tradeoff
#' The function `precision()` calculates a number of measures of the effect of shrinkage of the 
#' coefficients in relation to the “size” of the covariance matrix.

pridge <- precision(lridge) |> print()


plot(pridge, criteria = lridge$criteria)

plot(pridge, 
     labels = "df", label.prefix="df:")

plot(pridge, yvar="trace",
     labels = "df", label.prefix="df:")

#' ## Low-rank views
#' The `pca` method transforms a "ridge" object from parameter space to PCA space
#' defined by the SVD of X

plridge <- pca(lridge) |> print()

traceplot(plridge, X="df", 
          cex.lab = 1.2, lwd=2)

#' ## View the last two dimensions
plot(plridge, variables=5:6, 
     fill = TRUE, fill.alpha=0.15, cex.lab = 1.5)
text(plridge$coef[, 5:6], 
     label = lambdaf, 
     cex=1.5, pos=4, offset=.1)

#' ## Add variable vectors
biplot(plridge, radius=0.5, 
       ref=FALSE, asp=1, 
       var.cex=1.15, cex.lab=1.3, col=clr,
       fill=TRUE, fill.alpha=0.15, 
       prefix="Dimension ")
text(plridge$coef[,5:6], lambdaf, pos=2, cex=1.3)








