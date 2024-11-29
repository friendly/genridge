#' ---
#' title: Longley data example 
#' author: "Michael Friendly"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'      theme: journal
#'      code_download: true
#' ---

#+ echo=FALSE
knitr::opts_chunk$set(warning=FALSE, 
                      message=FALSE, 
                      R.options=list(digits=4))

#' ## Load packages
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
#' The ridge shrinkage factors are specified by `lambda` (aka "k")
lambda <- c(0, 0.002, 0.005, 0.01, 0.02, 0.04, 0.08)
lridge <- ridge(Employed ~ GNP + Unemployed + Armed.Forces + 
                  Population + Year + GNP.deflator, 
                data=longley, lambda=lambda)
print(lridge, digits = 2)

#' ## The `"ridge"` object
#' The main thing is the coefficients, (returned by `coef()`) but the `ridge` object contains
#' info used by other functions so it doesn't have to be computed again.
#'  
names(lridge)

#' ## Standard, univariate trace plot
#' `traceplot()` plots estimated coefficients against the ridge constant, λ
traceplot(lridge, 
          X = "lambda",
          xlab = "Ridge constant (λ)",
          xlim = c(-0.02, 0.08), 
          cex.lab=1.25)

#' ## Better: plot against effective degrees of freedom
#' Shrinkage can be better understood in terms of the equivalent number of degrees of freedom.
#' Note that these values are more nearly linearly spaced.
traceplot(lridge, 
          X = "df",
          xlim = c(4, 6.2), cex.lab=1.25)

#' ## Even BETTER: Bivariate ridge trace plot
#' Shows bias _and_ precision by confidence ellipses for the coefficients

clr <-  c("black", "red", "brown", "darkgreen","blue", "cyan4", "magenta")
pch <- c(15:18, 7, 9, 12)
lambdaf <- c(expression(~widehat(beta)^OLS), as.character(lambda[-1]))

#' Plot GNP vs. Unemployed
plot(lridge, variables=c(1,2),
     radius=0.5, cex.lab=1.5, col=clr, 
     labels=NULL, fill=TRUE, fill.alpha=0.2)
text(lridge$coef[1,1], lridge$coef[1,2], 
     expression(~widehat(beta)^OLS), cex=1.5, pos=4, offset=.1)
text(lridge$coef[-1,c(1,2)], lambdaf[-1], pos=3, cex=1.3)

#' Plot GNP vs. Population
plot(lridge, variables=c(1,4),
     radius=0.5, cex.lab=1.5, col=clr, 
     labels=NULL, fill=TRUE, fill.alpha=0.2)
text(lridge$coef[1,1], lridge$coef[1,4], 
     expression(~widehat(beta)^OLS), cex=1.5, pos=4, offset=.1)
text(lridge$coef[-1,c(1,4)], lambdaf[-1], pos=3, cex=1.3)

#' ## Scatterplot matrix
#' See all bivariate pairwise views of effects of shrinkage. Some of the paths are monotonic, while others are not.
#' 
pairs(lridge, radius=0.5, diag.cex = 2, 
      fill = TRUE, fill.alpha = 0.1)

#' ## Visualizing the bias-variance tradeoff
#' The function `precision()` calculates a number of measures of the effect of shrinkage of the 
#' coefficients in relation to the “size” of the covariance matrix.

pridge <- precision(lridge) |> print(digits = 3)

#' ## Plot method for `"precision"` objects
#' `genridge` recently gained a `plot()` method to plot a measures of shrinkage ("bias") against a measure of
#' variance of the covariance matrices (inverse "precision")
plot(pridge, criteria = lridge$criteria)

#' Label the points with `df` and compare `"det"` = log |Var(β)| vs. `"trace"` = tr(Var(β)). These have quite different trajectories.
#+ fig.show="hold", fig.width=9
op <- par(mar = c(4,4,1,1)+ .5, mfrow = c(1,2))
plot(pridge, 
     labels = "df", label.prefix="df:")

plot(pridge, yvar="trace",
     labels = "df", label.prefix="df:")
par(op)

#' ## Low-rank views
#' The `pca` method transforms a `"ridge"` object from parameter space to PCA space
#' defined by the SVD of X. It returns an object of class `"pcaridge"`.
#' Note there is very little shrinkage in the first 4 dimensions.

plridge <- pca(lridge) |> print(digits = 3)

traceplot(plridge, X="df", 
          cex.lab = 1.2, lwd=2)

#' ## View the last two dimensions
plot(plridge, variables=5:6, 
     fill = TRUE, fill.alpha=0.15, cex.lab = 1.5)
text(plridge$coef[, 5:6], 
     label = lambdaf, 
     cex=1.5, pos=4, offset=.1)

#' ## Add variable vectors: `biplot()` method
#' The `biplot()` method produces a 2D plot of the covariance ellipses for two dimensions in PCA space,
#' defaulting to the last two dimensions where things get interesting.
#' Variable vectors show the projection of the variables into this space, facilitating interpretation.
#' 
biplot(plridge, radius=0.5, 
       ref=FALSE, asp=1, 
       var.cex=1.15, cex.lab=1.3, col=clr,
       fill=TRUE, fill.alpha=0.15, 
       prefix="Dimension ")
text(plridge$coef[,5:6], lambdaf, pos=2, cex=1.3)








