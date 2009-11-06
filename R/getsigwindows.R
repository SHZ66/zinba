getsigwindows=function(file,formula,threshold=.01,peakconfidence=.8,priorpeakprop=.15,winout,coordout,tol=10^-5,getPeakRefine=1,method='pscl'){
	time.start <- Sys.time()
	library(qvalue)
        #library(pscl)
	library(MASS)
        options(scipen=999)
        
        #print(paste("Using",method,"method for detection",sep=" "))
        #print(paste("winout is",winout,"coordout is",coordout,sep=" "))
	###### USER INPUT############################
	if(method=='pscl'){
            files = unlist(strsplit(file,";"))
            for(i in 1:length(files)){
                data=read.table(files[i], header=TRUE)
                mf <- model.frame(formula=formula, data=data)
                X <- model.matrix(attr(mf, "terms"), data=mf)
                if(i == 1){
                        a=zeroinfl(formula, data=data,dist='negbin', EM=TRUE)
                }else{
                        a=zeroinfl(formula, data=data,dist='negbin', EM=TRUE,start=param)
                }
		
		q25=quantile(data$exp_count, 0.25)
                leverage=hat(X, intercept=FALSE)
                fdrlevel=threshold
                standardized=residuals(a)/sqrt(1-leverage)
                pval=1-pnorm(as.matrix(standardized))
                fdr=qvalue(pval)
                numpeaks=length(which(fdr[[3]]<fdrlevel))
                minresid=min(standardized[which(fdr[[3]]<fdrlevel)])
                sigpeaks=cbind(data[which(fdr[[3]]<fdrlevel),], fdr[[3]][which(fdr[[3]]<fdrlevel)], standardized[which(fdr[[3]]<fdrlevel)])
                colnames(sigpeaks)[c(dim(sigpeaks)[2]-1, dim(sigpeaks)[2])]=c('q-value', 'residual')
                param=list(count=a$coefficients$count, zero=a$coefficients$zero, theta=a$theta)

                lineBLANK=''
                line0=paste('For ',as.character(files[i]),sep='')
                line1='|Selection Summary|'
                line2=paste('Selected number of peaks: ', as.character(numpeaks), sep='')
                line3=paste('Minimum Standardized Residual Value of peaks: ', as.character(minresid), sep='')

    ### PRINT SIGNIFICANT WINDOWS
                print(c(lineBLANK,line0,line1,line2,line3,lineBLANK))
                if(file.exists(winout)){
                    write.table(sigpeaks,winout,quote=F,sep="\t",row.names=F,col.names=F,append=T)
                }else{
                    write.table(sigpeaks,winout,quote=F,sep="\t",row.names=F)
                }

    ### FORMAT PEAK COORDINATE DATA
                if(getPeakRefine == 1){
                    peakID=paste(sigpeaks$chromosome,sigpeaks$start,sigpeaks$stop,sep=":")
                    coordinates=cbind(peakID,as.character(sigpeaks$chromosome),sigpeaks$start,sigpeaks$stop,((sigpeaks$exp_count>q25)^2),"+")
                    write.table(coordinates,coordout,quote=F,append=TRUE,sep="\t",row.names=F,col.names=F)
                }
	    }
	}else if(method=='mixture'){
            files = unlist(strsplit(file,";"))
            for(i in 1:length(files)){
                data=read.table(files[i], header=TRUE)
		q25=quantile(data$exp_count, 0.25)
                mf <- model.frame(formula=formula, data=data)
                X <- model.matrix(attr(mf, "terms"), data=mf)
		
		logsumexp=function(v){
			if(any(is.infinite(v))){
				stop("infinite value in v\n")
			}
			if(length(v)==1){ return(v[1]) }
			sv  = sort(v, decreasing=TRUE)
			res = sum(exp(sv[-1] - sv[1]))
			lse = sv[1] + log(1+res)
			lse
		}
          
		loglikfun=function(parms){
			mu1=exp(X%*%parms[1:kx])
			mu2=exp(X%*%parms[(kx+1):(kx+kz)])
			theta1 <- exp(parms[(kx+kz+kz)+1])
		        theta2 <- exp(parms[(kx+kz+kz)+2])	
			prob0 <- as.vector(linkinv(Z %*% parms[(kx+kz+1):(kx+kz+kz)]))
			prop1=parms[kx+kz+kz+3]
			prop2=parms[kx+kz+kz+4]
			loglik0=log(prob0+prop1*dnbinom(0, size = theta1, mu = mu1)+prop2*dnbinom(0, size = theta2, mu = mu2))
			loglik1=log(prop1*dnbinom(Y, size = theta1, mu = mu1)+prop2*dnbinom(Y, size = theta2, mu = mu2))
			NAs=which(loglik1==-Inf | loglik1==-Inf)
			if(length(NAs>0)){
				loglik1[NAs]=apply(cbind(log(prop1)+dnbinom(Y[NAs], size = theta1, mu = mu1[NAs],log=TRUE),log(prop2)+dnbinom(Y[NAs], size = theta2, mu = mu2[NAs], log=TRUE)), 1, logsumexp)
			}
			loglik=sum(loglik0*Y0+loglik1*Y1)
			loglik
		}
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
			model_zero <- suppressWarnings(glm.fit(Z, as.integer(Y0),family = binomial(link = linkstr)))
			prop0=sum( model_zero$fitted)/n
	
			#starting params for count componenets
			if(i == 1){			
				prop2=priorpeakprop
			}
			prop1=1-prop0-prop2
			odY = order(Y)
			n1  = round(length(Y) * (1 - prop2))
			priorCOUNTweight=rep(1, length(Y))      
			priorCOUNTweight[odY[1:n1]]=0
	
			model_count1 <- suppressWarnings(glm.fit(X, Y, family = poisson(), weights = (1-priorCOUNTweight)))
			model_count2 <- suppressWarnings(glm.fit(X, Y, family = poisson(), weights = (priorCOUNTweight)))
	
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
			NAs=which(probi1=='NaN'| probi2=='NaN')			
			if(length(NAs>0)){
				probi1[NAs]=0
				probi2[NAs]=1
			}	
		ll_new <- loglikfun(c(start$count1, start$count2, start$zero, log(start$theta1), log(start$theta2), prop1, prop2))

		ll_old <- 2 * ll_new      

		if(!require("MASS")) {
			ll_old <- ll_new
			warning("EM estimation of starting values not available")
		}
		ll=matrix(0, 1, 10000)
		ll[1]=ll_new
		i=2
		while(abs((ll_old - ll_new)/ll_old) > tol) {
#		        print(ll_new)
			ll_old <- ll[max(1, i-10)]
			 prop1=sum(probi1)/n
			 prop2=sum(probi2)/n
			 #updated values for parameters of component means
		         
			model_zero <- suppressWarnings(glm.fit(Z, probi0, family = binomial(link = linkstr), start = start$zero))
			model_count1 <- glm.nb(Y ~ 0 + X, weights = (probi1),start = start$count1, init.theta = start$theta1)
			model_count2 <- glm.nb(Y ~ 0 + X, weights = (probi2),start = start$count2, init.theta = start$theta2)
			start <- list(count1 = model_count1$coefficients, count2 = model_count2$coefficients,theta1 = model_count1$theta, theta2 = model_count2$theta, zero = model_zero$coefficients)

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

			ll_new <- loglikfun(c(start$count1, start$count2, start$zero, log(start$theta1), log(start$theta2), prop1, prop2))

			ll[i]=ll_new
			i=i+1 
		}
		numpeaks=length(which(probi2>peakconfidence))	

		minresid='NA'
		sigpeaks=cbind(data[which(probi2>peakconfidence),],probi2[probi2>peakconfidence])
		colnames(sigpeaks)[dim(sigpeaks)[2]]='peakprob'

                lineBLANK = ''
                line0 = paste("Processing ",files[i])
		line1='|Selection Summary|'
		line2=paste('Selected number of peaks: ', as.character(numpeaks), sep='')
		line3=paste('Minimum Standardized Residual Value of peaks: ', as.character(minresid), sep='')
        ### PRINT SIGNIFICANT WINDOWS
		print(c(lineBLANK,line0,line1,line2,line3,lineBLANK))
                if(file.exists(winout)){
                    write.table(sigpeaks,winout,quote=F,sep="\t",row.names=F,col.names=F,append=T)
                }else{
                    write.table(sigpeaks,winout,quote=F,sep="\t",row.names=F)
                }

    ### FORMAT PEAK COORDINATE DATA
                if(getPeakRefine == 1){
                    peakID=paste(sigpeaks$chromosome,sigpeaks$start,sigpeaks$stop,sep=":")
                    coordinates=cbind(peakID,as.character(sigpeaks$chromosome),sigpeaks$start,sigpeaks$stop,((sigpeaks$exp_count>q25)^2),"+")
                    write.table(coordinates,coordout,quote=F,append=TRUE,sep="\t",row.names=F,col.names=F)
                }
            }
	}
	time.end <- Sys.time()
	print(difftime(time.end,time.start))
}
