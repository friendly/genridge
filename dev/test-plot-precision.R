# test plots for precision

# added one more lambda
lambda <- c(0, 0.001, 0.005, 0.01, 0.02, 0.04, 0.08)
lambdaf <- c(expression(~widehat(beta)^OLS), lambda[-1])
lridge <- ridge(Employed ~ GNP + Unemployed + Armed.Forces + 
                  Population + Year + GNP.deflator, 
                data=longley, lambda=lambda)
pdat <- precision(lridge)

clr <-  c("black", "red", "darkgreen","blue", "cyan4", "magenta", "brown")
#pch <- c(15:18, 7, 9)

#op <- par(mar=c(4, 4, 1, 1) + 0.2)

library(splines)
with(pdat, {
  plot(norm.beta, det, type="b", 
       cex.lab=1.25, pch=16, cex=1.5, col=clr, lwd=2,
       xlab='shrinkage: ||b|| / max(||b||)',
       ylab='variance: log |Var(b)|')
  text(norm.beta, det, 
       labels = lambdaf, 
       cex=1.25, pos=c(rep(2,length(lambda)-1),4), xpd = TRUE)
  text(min(norm.beta), max(det), 
       labels = "log |Variance| vs. Shrinkage", 
       cex=1.5, pos=4)
})
# find locations for optimal shrinkage criteria
mod <- lm(cbind(det, norm.beta) ~ bs(lambda, df=5), 
          data=pdat)
x <- data.frame(lambda=c(lridge$kHKB, 
                         lridge$kLW))
fit <- predict(mod, x)
points(fit[,2:1], pch=15, col=gray(.50), cex=1.6)
text(fit[,2:1], c("HKB", "LW"), pos=3, cex=1.5, col=gray(.50))

# --- trace -----

with(pdat, {
  plot(norm.beta, trace, type="b", 
       cex.lab=1.25, pch=16, cex=1.5, col=clr, lwd=2,
       xlab='shrinkage: ||b|| / max(||b||)',
       ylab='variance: trace[Var(b)]')
  text(norm.beta, trace, 
       labels = lambdaf, 
       cex=1.25, pos=c(rep(2,length(lambda)-1),4), xpd = TRUE)
  text(min(norm.beta), max(trace), 
       labels = "trace(Variance) vs. Shrinkage", 
       cex=1.5, pos=4)
})
# find locations for optimal shrinkage criteria
mod <- lm(cbind(trace, norm.beta) ~ bs(lambda, df=5), 
          data=pdat)
x <- data.frame(lambda=c(lridge$kHKB, 
                         lridge$kLW))
fit <- predict(mod, x)
points(fit[,2:1], pch=15, col=gray(.50), cex=1.6)
text(fit[,2:1], c("HKB", "LW"), pos=3, cex=1.5, col=gray(.50))

# --- max.eig --------

with(pdat, {
  plot(norm.beta, max.eig, type="b", 
       cex.lab=1.25, pch=16, cex=1.5, col=clr, lwd=2,
       xlab='shrinkage: ||b|| / max(||b||)',
       ylab='variance: max eigenval Var(b)')
  text(norm.beta, max.eig, 
       labels = lambdaf, 
       cex=1.25, pos=c(rep(2,length(lambda)-1),4))
  text(min(norm.beta), max(max.eig), 
       labels = "max eigval(Variance) vs. Shrinkage", 
       cex=1.5, pos=4)
})
# find locations for optimal shrinkage criteria
mod <- lm(cbind(max.eig, norm.beta) ~ bs(lambda, df=5), 
          data=pdat)
x <- data.frame(lambda=c(lridge$kHKB, 
                         lridge$kLW))
fit <- predict(mod, x)
points(fit[,2:1], pch=15, col=gray(.50), cex=1.6)
text(fit[,2:1], c("HKB", "LW"), pos=3, cex=1.5, col=gray(.50))

