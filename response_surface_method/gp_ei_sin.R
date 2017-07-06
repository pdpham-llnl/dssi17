library(laGP)                                      ## for GP fitting
library(plgp, quietly=TRUE, warn.conflicts=FALSE)  ## for EI calculations
library(mvtnorm)                                   ## for visualization

## testing initialization code
start <- 6 ## (don't go less than 6)
X <- matrix(seq(0,2*pi,length=start), ncol=1)
Z <- sin(X)

## deliberately small d=1 value to make for a better illustration
gpi <- newGP(X, Z, 1, 0.00001, dK=TRUE)
## add mleGP after updateGP below to speed convergenves

## testing prediction code
XX <- matrix(seq(-1,2*pi+1, length=499), ncol=ncol(X))

end <- 20
for(t in (start+1):end) {

  ## calculate EI
  outp <- predGP(gpi, XX)  ## lite = TRUE is faster but can't use rmvt below
  tmat <- data.frame(m=outp$mean, sd=sqrt(diag(outp$Sigma)), df=outp$df)
  eis <- calc.eis(tmat, min(Z))

  ## calculate next point
  m <- which.max(eis)
  Xnew <- matrix(XX[m,], ncol=ncol(X))
  Znew <- sin(Xnew)

  ## plot progress
  par(mfrow=c(2,1), mar=c(4,4,0.5,1))
  N <- 100
  ZZ <- rmvt(N, outp$Sigma, outp$df) + t(matrix(rep(outp$mean, N), ncol=N))
  matplot(XX, t(ZZ), col="gray", lwd=0.5, lty=1, type="l",
          ylab="Y", xlab="", ylim=c(-1.3,1.5))
  points(X, Z, pch=19, cex=1.5)
  lines(XX[,1], outp$mean, lwd=2)
  sd <- sqrt(diag(outp$Sigma))
  lines(XX[,1], outp$mean + qt(0.95, outp$df)*sd, col=2, lty=2, lwd=2)
  lines(XX[,1], outp$mean + qt(0.05, outp$df)*sd, col=2, lty=2, lwd=2)
  points(Xnew, Znew, col=3, pch=18, cex=2)
  plot(XX, eis, type="l", ylab="EI", xlab="X", ylim=c(0,0.25), lwd=2, col="blue")

  ## update GP fit
  readline("press RETURN to continue: ")
  X <- rbind(X, Xnew)
  Z <- c(Z, Znew)
  XX <- matrix(XX[-m,], ncol=ncol(X))
  updateGP(gpi, Xnew, Znew, verb=1)
  ## mleGP(gpi, tmax=10)
}

## clean up
deleteGP(gpi)
