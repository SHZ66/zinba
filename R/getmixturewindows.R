getmixturewindows=function(formula, data, outputpath,priorpeakprop=.15, peakconfidence=.8){
time.start <- Sys.time()

#log-likelihood
loglikfun=function(parms){
	mu1=exp(X%*%parms[1:kx])
	mu2=exp(X%*%parms[(kx+1):(kx+kz)])
	theta1 <- exp(parms[(kx+kz+kz)+1])
        theta2 <- exp(parms[(kx+kz+kz)+2])	
	prob0 <- as.vector(linkinv(Z %*% parms[(kx+kz+1):(kx+kz+kz)]))
	prop1=parms[kx+kz+kz+3]
	prop2=parms[kx+kz+kz+4]
	loglik=sum(log(prob0*Y0+prop1*dnbinom(Y, size = theta1, mu = mu1)+prop2*dnbinom(Y, size = theta2, mu = mu2)))
	loglik
}


#data setup
mf <- model.frame(formula=formula, data=data)
X <- model.matrix(attr(mf, "terms"), data=mf)
Y <- model.response(mf)
Z=X
n <- length(Y)
kx <- NCOL(X)
kz <- NCOL(Z)
Y0 <- Y <= 0
Y1 <- Y > 0

linkstr <- 'logit'
linkobj <- make.link(linkstr)
linkinv <- linkobj$linkinv


#starting params for ze0 component
model_zero <- glm.fit(Z, as.integer(Y0),family = binomial(link = linkstr))
prop0=sum( model_zero$fitted)/n

#starting params for count componenets
prop2=priorpeakprop
prop1=1-prop0-prop2
odY = order(Y)
n1  = round(length(Y) * (1 - prop2))
priorCOUNTweight=rep(1, length(Y))      
priorCOUNTweight[odY[1:n1]]=0

model_count1 <- glm.fit(X, Y, family = poisson(), weights = (1-priorCOUNTweight))
model_count2 <- glm.fit(X, Y, family = poisson(), weights = (priorCOUNTweight))

#start parameter vector
start <- list(count1 = model_count1$coefficients, count2 = model_count2$coefficients,zero = model_zero$coefficients)
start$theta1 <- 1
start$theta2 <- 1

#starting prior probs
mui1  <- model_count1$fitted
mui2  <- model_count2$fitted
probi0 <- model_zero$fitted

probi0=probi0/(probi0+prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ prop2*dnbinom(Y, size = start$theta2, mu = mui2))
probi0[Y1]=0
probi1  <- prop1*dnbinom(Y, size = start$theta1, mu = mui1)/(probi0*Y0+prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ prop2*dnbinom(Y, size = start$theta2, mu = mui2))
probi2  <- prop2*dnbinom(Y, size = start$theta2, mu = mui2)/(probi0*Y0+prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ prop2*dnbinom(Y, size = start$theta2, mu = mui2))

ll_new <- loglikfun(c(start$count1, start$count2, start$zero, log(start$theta1), log(start$theta2), prop1, prop2))

ll_old <- 2 * ll_new      

if(!require("MASS")) {
	ll_old <- ll_new
	warning("EM estimation of starting values not available")
}
ll=matrix(0, 1, 10000)
ll[1]=ll_new
i=2
while(abs((ll_old - ll_new)/ll_old) > 10^-5) {
	print(ll_new)
	ll_old <- ll[max(1, i-10)]
	 prop1=sum(probi1)/n
	 prop2=sum(probi2)/n
	 #updated values for parameters of component means
         
	model_zero <- glm.fit(Z, probi0, family = binomial(link = linkstr), start = start$zero)
	model_count1 <- glm.nb(Y ~ 0 + X, weights = (probi1),start = start$count1, init.theta = start$theta1)
	model_count2 <- glm.nb(Y ~ 0 + X, weights = (probi2),start = start$count1, init.theta = start$theta1)
	start <- list(count1 = model_count1$coefficients, count2 = model_count2$coefficients,theta1 = model_count1$theta, theta2 = model_count2$theta, zero = model_zero$coefficients)

	mui1  <- model_count1$fitted
	mui2  <- model_count2$fitted
	probi0 <- model_zero$fitted

	probi0=probi0/(probi0+prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ prop2*dnbinom(Y, size = start$theta2, mu = mui2))
	probi0[Y1]=0
	probi1  <- prop1*dnbinom(Y, size = start$theta1, mu = mui1)/(probi0*Y0+prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ prop2*dnbinom(Y, size = start$theta2, mu = mui2))
	probi2  <- prop2*dnbinom(Y, size = start$theta2, mu = mui2)/(probi0*Y0+prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ prop2*dnbinom(Y, size = start$theta2, mu = mui2))
	ll_new <- loglikfun(c(start$count1, start$count2, start$zero, log(start$theta1), log(start$theta2), prop1, prop2))
	ll[i]=ll_new
	i=i+1 
}
time.end <- Sys.time()
print(time.start-time.end)
print(start)
peaks=cbind(data[which(probi2>peakconfidence),],probi2[probi2>peakconfidence])

write.table(peaks,outputpath)
}
