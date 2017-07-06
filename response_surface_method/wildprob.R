set.seed(1)
graphics.off()

## aimprob:
##
## evaulate 2-d test problem returning objective and constraints

aimprob <- function(x)
{
  f <- sum(x)
  c1 <- 1.5 - x[1] - 2*x[2] - 0.5*sin(2*pi*(x[1]^2 - 2*x[2]))
  c2 <- x[1]^2 + x[2]^2-1.5
  return(list(obj=f, c=c(c1,c2)))
}

## aimprob.AL:
##
## augmented Lagrangian version

aimprob.AL <- function(x, B, lambda, rho)
{
	if(any(x < B[,1]) | any(x > B[,2])) return(Inf)
	fc <- aimprob(x)
	al <- fc$obj + lambda %*% fc$c + rep(1/(2*rho), 2) %*% pmax(0, fc$c)^2
	evals <<- evals + 1
	return(al) 
}

## bounding box
B <- matrix(c(rep(0,2),rep(1,2)),ncol=2)

## multiple restarts
while(1) {

## keep track of number of aimprob evaluations
evals <- 0

## initial AL values
lambda <- rep(0, ncol(B))
rho <- 1/2

## set up a plot
plot(0, 0, type="n", xlim=B[1,], ylim=B[2,], xlab="x1", ylab="x2")
x <- seq(0,1, length=100)
X <- expand.grid(x, x)
c1 <- 1.5 - X[,1] - 2*X[,2] - 0.5*sin(2*pi*(X[,1]^2 - 2*X[,2]))
contour(x, x, matrix(c1, ncol=length(x)), nlevels=1, levels=0, 
	drawlabels=FALSE, add=TRUE, lwd=2)
c2 <- X[,1]^2 + X[,2]^2 - 1.5
contour(x, x, matrix(c2, ncol=length(x)), nlevels=1, levels=0, 
	drawlabels=FALSE, add=TRUE, lwd=2)

## random starting point
start <- runif(2)
fc <- aimprob(start)

## add starting point to plot
points(start[1], start[2], cex=0.75, pch=19,
	col={if(all(fc$c <= 0)) "black" else "red"})	
readline("press RETURN to continue: ")

## for keeping track of best
best.obj <- Inf
best <- NULL

## iterations of Augmented Lagrangian
for(i in 1:10) {

	## solve sub-problem
	out <- optim(start, aimprob.AL, control=list(maxit=15),
		B=B, lambda=lambda, rho=rho)

	## extract x^* that was found, and plot progress
	start <- out$par
	fc <- aimprob(start)
	points(start[1], start[2], cex=0.75, 
		col={if(all(fc$c <= 0)) "black" else "red"})
	readline("press RETURN to continue: ")

	## check if new best valid
	if(all(fc$c <= 0) && fc$obj < best.obj) {
		best.obj <- fc$obj
		best <- start
	}

	## update augmented Lagrangian
	lambda <- pmax(0, lambda + (1/rho)*fc$c)
	if(any(fc$c > 0)) rho = rho/2
}

## put cross-hairs on the best
abline(v=best[1], col="green", lty=2)
abline(h=best[2], col="green", lty=2)
## total number of evaluations
title(paste(evals, "aimprob evaluations"))

## repeat?
readline("demo done; press RETURN to repeat: ")
}