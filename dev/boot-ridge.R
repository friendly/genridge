#################################
## Bootstrapping ridge regression
#################################

if (!require(genridge)) {install.packages("genridge"); library(genridge)}
library(car)

longley.y <- longley[, "Employed"]
longley.X <- data.matrix(longley[, c(2:6,1)])

lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
lambdaf <- c("", ".005", ".01", ".02", ".04", ".08")
lridge <- ridge(longley.y, longley.X, lambda=lambda)
col <- c("black", rainbow(5, start=.6, end=.1))

# function to bootstrap ridge regression
# NB:  need to trap cases which give a singular result
boot.ridge <- function(data, indices, lambda) {
	data <- data[indices,]
	y <- data[,1]
	x <- data[,-1]
	mod <- try(ridge(y, x, lambda=lambda, svd=FALSE))
	if (!inherits(mod, "try-error")) coef(mod) else NA
}


set.seed(12313147)
library(boot)
Longley <- data.frame(Employed=longley.y, longley.X)
B <- 800

lam <- lambda[c(1,3,6)]
clr <- col[c(1,3,6)]
lboot <- boot(Longley, boot.ridge, B, lambda=lam)

colnames(lboot$t) <- c(t(outer(colnames(lboot$t0), rownames(lboot$t0), paste, sep=':')))

# extract bootstrap estimates
# transform to a data.frame with separate column for lambda
nl <- length(lam)
cols <- rep(1:6, nl) 
lbdata <- matrix(0,nrow=nl*nrow(lboot$t), ncol=6)
for (i in 1:nrow(lboot$t)) {
	lbdata[(nl*i-2):(nl*i),] <- matrix(lboot$t[i,], nrow=3)
}
lbdata<- as.data.frame(lbdata)
colnames(lbdata) <- colnames(lboot$t0)
lbdata$lambda <- rep(lam, length=nrow(lbdata))
str(lbdata)


lbt <- lboot$t
wh <- c(1,4)
xlab <- colnames(lboot$t0)[1]
ylab <- colnames(lboot$t0)[2]
cm0 <- colMeans(lbt[,wh])
cm1 <- colMeans(lbt[,1+wh])
cm2 <- colMeans(lbt[,2+wh])

# levels=.12 -> radius=.50 
#(radius <- sqrt(2* qf(.12, 2, 799)))

##############################
# Figure 10a
op <- par(mar=c(4, 4, 1, 1) + 0.4)
dataEllipse(lbt[,wh], add=FALSE, levels=0.12, col=clr[1], 
	plot.points=TRUE, cex=0.4, lwd=3,
	xlim=c(-8,4), ylim=c(-2.5,-0.5),
	xlab="GNP", ylab="Unemployed", cex.lab=1.25
	)
dataEllipse(lbt[,1+wh], plot.points=FALSE, 
	lwd=3,
	add=TRUE, levels=0.12, col=clr[2])
dataEllipse(lbt[,2+wh], 
	plot.points=FALSE, 
#	plot.points=TRUE, cex=0.3, lwd=3,
  add=TRUE, levels=0.12, col=clr[3])
lines(rbind(cm0, cm1, cm2))

text(cm0[1], cm0[2]-.02, expression(~widehat(beta)^OLS), pos=3, cex=1.35)
text(cm1[1], cm1[2]+.08, lam[2], pos=3, cex=1.25)
text(cm2[1], cm2[2]+.05, lam[3], pos=3, cex=1.25)
text(-8, -0.6, "Bootstrap RR: data ellipses", cex=1.5, pos=4)
par(op)

dev.copy2eps(file="ridge-boot1.eps")
dev.copy2pdf(file="ridge-boot1.pdf")


library(KernSmooth)
dest0 <- bkde2D(lbt[,wh], bandwidth=.4, gridsize=c(81,81), range.x=list(c(-16, 8), c(-5,5)))
dest1 <- bkde2D(lbt[,1+wh], bandwidth=.4, gridsize=c(81,81), range.x=list(c(-16, 8), c(-5,5)))
dest2 <- bkde2D(lbt[,2+wh], bandwidth=.4, gridsize=c(81,81), range.x=list(c(-16, 8), c(-5,5)))
op <- par(mar=c(4, 4, 1, 1) + 0.4)
contour(dest0$x1, dest0$x2, dest0$fhat, nlevels=4,
	xlim=c(-8,4), ylim=c(-3,0), cex.lab=1.25,
	xlab="GNP", ylab="Unemployed" 	)
contour(dest1$x1, dest1$x2, dest1$fhat, nlevels=4, add=TRUE, col=clr[2])
contour(dest2$x1, dest2$x2, dest2$fhat, nlevels=4, add=TRUE, col=clr[3])
lines(rbind(cm0, cm1, cm2))
points(rbind(cm0, cm1, cm2), col=clr, pch=16, cex=1.5)
text(cm0[1], cm0[2], expression(~widehat(beta)^OLS), pos=3, cex=1.25)
text(cm1[1], cm1[2], lam[2], pos=3, cex=1.25)
text(cm2[1], cm2[2], lam[3], pos=3, cex=1.25)
text(-8, -0.1, "Bootstrap RR: 2D Kernel density", cex=1.5, pos=4)
par(op)
dev.copy2eps(file="ridge-boot2.eps")
dev.copy2pdf(file="ridge-boot2.pdf")

# use my modification contourf for filled contours
source("contourf.R")

##############################
# Figure 10b
op <- par(mar=c(4, 4, 1, 1) + 0.4)
colramp = colorRampPalette(c("white", clr[1]))
cr <- paste(colramp(5), "80", sep="") # make transparent
contourf(dest0$x1, dest0$x2, dest0$fhat, nlevels=4,
	xlim=c(-8,4), ylim=c(-3,0), cex.lab=1.25,
	xlab="GNP", ylab="Unemployed", fill.col=cr)

colramp = colorRampPalette(c("white", clr[2]))
cr <- paste(colramp(5), "A0", sep="") # make transparent
contourf(dest1$x1, dest1$x2, dest1$fhat, nlevels=4, add=TRUE, col=clr[2], fill.col=cr)

colramp = colorRampPalette(c("white", clr[3]))
cr <- paste(colramp(5), "A0", sep="") # make transparent
contourf(dest2$x1, dest2$x2, dest2$fhat, nlevels=4, add=TRUE, col=clr[3], fill.col=cr)
lines(rbind(cm0, cm1, cm2))
points(rbind(cm0, cm1, cm2), col=clr, pch=16, cex=1.5)
text(cm0[1], cm0[2], expression(~widehat(beta)^OLS), pos=3, cex=1.25)
text(cm1[1], cm1[2], lam[2], pos=3, cex=1.25)
text(cm2[1], cm2[2], lam[3], pos=3, cex=1.25)
text(-8, -0.1, "Bootstrap RR: 2D Kernel density", cex=1.5, pos=4)
par(op)

#dev.copy2eps(file="ridge-boot2f.eps")
dev.copy2pdf(file="ridge-boot2f.pdf")

