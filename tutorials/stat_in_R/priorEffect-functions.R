###########################
#author: Ana Paula Sales  #
#last edited: 07/06/2016  #
###########################


############################
# Prior Effect
#TRUE distribution: x ~ N(mu,sd)
#PRIOR: mu ~ N(mu0,sd0)
############################
sim.conjNormKnownVar = function(prior.mean = 3, 
                                prior.sd=1, 
                                true.mean=0, 
                                known.sd=1, 
                                max.y = .7, 
                                n.samples =1, 
                                minIterationToPlot =0,
                                numIterationsToSkip =0,
                                interactive = FALSE, #if TRUE, waits for keyboard input before each new plot
                                x = NULL,  #if NULL, will sample x from N(true.mean, true.sd) 
                                save.plots = FALSE, 
                                output.dir='./', 
                                draw.true.mean.line=TRUE){
  # x range of plot
  max.x = max(prior.mean, true.mean) + 3*max(prior.sd, known.sd)
  min.x = min(prior.mean, true.mean) - 3*max(prior.sd, known.sd)
  plot.sup = seq(min.x,max.x,length=1000)
  
  #prior density 
  prior.dens = dnorm(plot.sup,prior.mean,prior.sd)
  
  #if data not provided, sample data from true distribution
  if(is.null(x)){
    x = rnorm(n.samples, true.mean, known.sd)
  }
   
  #if writing plots to file, check if directory exists and if not, create it
  if(save.plots){
    if(!file.exists(output.dir))
      dir.create(output.dir)
  }
  
  #iterate
  skipped = numIterationsToSkip
  for(i in 1:n.samples){
    post.mean = ( prior.mean/prior.sd^2 + sum(x[1:i])/known.sd^2) / (1/prior.sd^2 + i/known.sd^2)
    post.sd = sqrt( 1/(1/prior.sd^2 + i/known.sd^2))
    post.dens = dnorm(plot.sup, post.mean,post.sd)
    mle.dens = dnorm(plot.sup, mean(x[1:i]), known.sd)
    if(i > minIterationToPlot & skipped == numIterationsToSkip){
      skipped=1
      if(save.plots){
        pdf(paste(output.dir, sprintf("%05d.pdf",i),sep=''), width=10, height=5)
      }
      par(mfrow=c(1,2))
      #plot f(mu) ~ mu and f(mu|x) ~ mu
      plot(plot.sup, prior.dens, type="l", ylim = c(0,max.y),  xlab="Mu",ylab=" ", bty='n')#yaxt='n',  #plot prior
      lines(plot.sup, post.dens,type="l", col="green", lwd=2) #plot posterior
      if(draw.true.mean.line)
        abline(v=true.mean,col="red")
      legend("topright", c("Prior", paste0("Posterior N(", round(post.mean,d=2), ",", round(post.sd^2,d=2) ,")"), "True mean"), col=c("black","green", "red"), lty=1,bty='n')
      
      #plot f(x) ~ x 
      hist(x[1:i], freq=FALSE,   ylab=" ", xlab="X", xlim=range(plot.sup),ylim = c(0,max.y), border="purple",bty='n',main=paste0("n = ",i))#, breaks=8)
      points(x[1:i],rep(0,i),pch='x',col="purple")
      points(x[i],0,pch='x',col="black",cex=1.4)
      lines(plot.sup, mle.dens,col="blue",lwd=2) #map or mle dens
      legend("topright", c("Likelihood (MLE)", "Earlier data", "New obs"), col=c("blue","purple","black"), lty=c(1,0,0),pch=c(NA,'x','x'),bty='n')
      
      if(save.plots)
        dev.off()
      
      if(interactive)
        readline()
    }else{
      skipped = skipped+1
    }
  }
  #return(max(post.dens))
}

############################
# Prior Effect
#beta-binomial
#beta, alpha >0, beta >0
############################
sim.betaBinom = function(alpha = 1, beta = 1, n.trials =10, n.successes=5){
  plot.sup = seq(0,1,length=1000)
  
  prior.dens = dbeta(plot.sup, alpha, beta)
  
  post.alpha = alpha + n.successes
  beta.post = beta + n.trials - n.successes
  
  plot(plot.sup, prior.dens, type="l", ylim = c(0,max.y), main = paste("Post.mean=", post.mean, " Post.sd=", post.sd), ylab="Density", xlab="mu") #plot prior
  
}



