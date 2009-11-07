startenrichment=function(range){
probs=seq(range[1], range[2],,5)
result=rep(0, length(probs))
a=Sys.time()
for(k in 1:length(probs)){
			model_zero <-.C("pglm_fit", family=as.integer(1), N=as.integer(length(Y)), M=as.integer(ncol(X)), y=as.double(Y0), prior=as.double(rep(1,n)), offset=as.double(rep(0,n)), X=as.double(unlist(X)),  stratum=as.integer(rep(1,n)),init=as.integer(0), rank=integer(1), Xb=double(n*ncol(X)), fitted=as.double((rep(1,n) * Y0 + 0.5)/(rep(1,n) + 1)), resid=double(n), weights=double(n),scale=double(1), df_resid=integer(1), theta=as.double(-1), package='zinba')  
			prop0=sum( model_zero$fitted)/n
	
			#starting params for count componenets
			prop2=probs[k]
			prop1=1-prop0-prop2
			odY = order(Y)
			n1  = round(length(Y) * (1 - prop2))
			priorCOUNTweight=rep(1-10^-10, length(Y))      
			priorCOUNTweight[odY[1:n1]]=10^-10
	
			model_count1 <- .C("pglm_fit", family=as.integer(2), N=as.integer(length(Y)), M=as.integer(ncol(XNB)), y=as.double(Y), prior=as.double(1-priorCOUNTweight), offset=as.double(rep(0,length(Y))), X=as.double(unlist(XNB)),  stratum=as.integer(rep(1,length(Y))),init=as.integer(1), rank=integer(1), Xb=double(length(Y)*ncol(XNB)), fitted=as.double(Y+(Y==0)/6), resid=double(length(Y)), weights=double(length(Y)),scale=double(1), df_resid=integer(1), theta=as.double(-1), package='zinba')  

			model_count2 <- .C("pglm_fit", family=as.integer(2), N=as.integer(length(Y)), M=as.integer(ncol(XNB)), y=as.double(Y), prior=as.double(priorCOUNTweight), offset=as.double(rep(0,length(Y))), X=as.double(unlist(XNB)),  stratum=as.integer(rep(1,length(Y))),init=as.integer(1), rank=integer(1), Xb=double(length(Y)*ncol(XNB)), fitted=as.double(Y+(Y==0)/6), resid=double(length(Y)), weights=double(length(Y)),scale=double(1), df_resid=integer(1), theta=as.double(-1), package='zinba')  

			#starting prior probs
			mui1  <- model_count1$fitted
			mui2  <- model_count2$fitted
			probi0 <- model_zero$fitted

		        #start mean vector
			start <- list(count1 = model_count1$fitted, count2 = model_count2$fitted,zero = model_zero$fitted, zerocoef=model_zero$coefficients)
			start$theta1 <- 1
			start$theta2 <- 1
	
	
			probi0=probi0/(probi0+prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ prop2*dnbinom(Y, size = start$theta2, mu = mui2))
			probi0[Y1]=0
			probi1  <- prop1*dnbinom(Y, size = start$theta1, mu = mui1)/(probi0*Y0+prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ prop2*dnbinom(Y, size = start$theta2, mu = mui2))
			probi2  <- prop2*dnbinom(Y, size = start$theta2, mu = mui2)/(probi0*Y0+prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ prop2*dnbinom(Y, size = start$theta2, mu = mui2))	
			NAs=which(probi1=='NaN'| probi2=='NaN')			
			if(length(NAs>0)){
				probi1[NAs]=0
				probi2[NAs]=1
			}	
		ll_new <- loglikfun(list(start=start, prop1=prop1, prop2=prop2))

		ll_old <- 2 * ll_new      

		if(!require("MASS")) {
			ll_old <- ll_new
			warning("EM estimation of starting values not available")
		}
		ll=matrix(0, 1, 10000)
		ll[1]=ll_new
		i=2
		while(i<4) {
#		        print(ll_new)
			ll_old <- ll[max(1, i-10)]
			 prop1=sum(probi1)/n
			 prop2=sum(probi2)/n
			 #updated values for parameters of component means
		         
			model_zero <- .C("pglm_fit", family=as.integer(1), N=as.integer(length(Y)), M=as.integer(ncol(X)), y=as.double(probi0), prior=as.double(rep(1,n)), offset=as.double(rep(0,n)), X=as.double(unlist(X)),  stratum=as.integer(rep(1,n)),init=as.integer(0), rank=integer(1), Xb=double(n*ncol(X)), fitted=as.double((rep(1,n) * probi0 + 0.5)/(rep(1,n) + 1)), resid=double(n), weights=double(n),scale=double(1), df_resid=integer(1), theta=as.double(-1), package='zinba')  

			model_count1 <- .C("pglm_fit", family=as.integer(0), N=as.integer(length(Y)), M=as.integer(ncol(XNB)), y=as.double(Y), prior=as.double(probi1), offset=as.double(rep(0,length(Y))), X=as.double(unlist(XNB)),  stratum=as.integer(rep(1,length(Y))),init=as.integer(1), rank=integer(1), Xb=double(length(Y)*ncol(XNB)), fitted=as.double(start$count1), resid=double(length(Y)), weights=double(length(Y)),scale=double(1), df_resid=integer(1), theta=as.double(start$theta1), package='zinba')  

			model_count2 <- .C("pglm_fit", family=as.integer(0), N=as.integer(length(Y)), M=as.integer(ncol(XNB)), y=as.double(Y), prior=as.double(probi2), offset=as.double(rep(0,length(Y))), X=as.double(unlist(XNB)),  stratum=as.integer(rep(1,length(Y))),init=as.integer(1), rank=integer(1), Xb=double(length(Y)*ncol(XNB)), fitted=as.double(start$count2), resid=double(length(Y)), weights=double(length(Y)),scale=double(1), df_resid=integer(1), theta=as.double(start$theta2), package='zinba')  
 

			start <- list(count1 = model_count1$fitted, count2 = model_count2$fitted,zero = model_zero$fitted, zerocoef=model_zero$coefficients)
			start$theta1 <- model_count1$theta
			start$theta2 <- model_count2$theta

			mui1  <- model_count1$fitted
			mui2  <- model_count2$fitted
			probi0 <- model_zero$fitted

			probi0=probi0/(probi0+prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ prop2*dnbinom(Y, size = start$theta2, mu = mui2))
			probi0[Y1]=0
			probi1  <- prop1*dnbinom(Y, size = start$theta1, mu = mui1)/(probi0*Y0+prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ prop2*dnbinom(Y, size = start$theta2, mu = mui2))
			probi2  <- prop2*dnbinom(Y, size = start$theta2, mu = mui2)/(probi0*Y0+prop1*dnbinom(Y, size = start$theta1, mu = mui1)+ prop2*dnbinom(Y, size = start$theta2, mu = mui2))
			NAs=which(probi1=='NaN'| probi2=='NaN')			
			if(length(NAs>0)){
				probi1[NAs]=0
				probi2[NAs]=1
			}

			ll_new <- loglikfun(list(start=start, prop1=prop1, prop2=prop2))

			ll[i]=ll_new
			i=i+1 
		}
result[k]=ll_new
}

return(probs[which.max(result)])
}
