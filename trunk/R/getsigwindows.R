getsigwindows=function(file,formula,formulaE,threshold=.01,peakconfidence=.8,winout,printFullOut=0,tol=10^-5,method='pscl',initmethod, diff=0){
    time.start <- Sys.time()
    suppressPackageStartupMessages(library(qvalue))
    suppressPackageStartupMessages(library(quantreg))
    library(MASS)
    options(scipen=999)

    winfile=NULL
    if(method=='pscl'){
        files = unlist(strsplit(file,";"))
        for(i in 1:length(files)){
            data=read.table(files[i], header=TRUE)
            chrm = data$chromosome[1]
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
            standardized=residuals.zeroinfl(a)/sqrt(1-leverage)
            pval=1-pnorm(as.matrix(standardized))
            fdr=qvalue(pval)
            numpeaks=length(which(fdr[[3]]<fdrlevel))
	    if(printFullOut == 1){
		data=cbind(data, ((data$exp_count>q25)^2), fdr[[3]], standardized)
		colnames(data)[c(dim(data)[2]-2, dim(data)[2]-1, dim(data)[2])]=c('q25','qvalue', 'residual')
	    }else{
		data=cbind(data[1:3],((data$exp_count>q25)^2), fdr[[3]], standardized)
		colnames(data)=c('chromosome','start','stop','q25','qvalue', 'residual')
	    }
            param=list(count=a$coefficients$count, zero=a$coefficients$zero, theta=a$theta)

### PRINT SIGNIFICANT WINDOWS
            print(paste('For ',files[i],', found ',as.character(numpeaks),' significant wins',sep=''))
            winfile = paste(winout,"_",chrm,".wins",sep="")
            if(file.exists(winfile)){
                write.table(data,winfile,quote=F,sep="\t",row.names=F,col.names=F,append=T)
            }else{
                write.table(data,winfile,quote=F,sep="\t",row.names=F)
            }
        }
        time.end <- Sys.time()
        print(difftime(time.end,time.start))
        return(winfile)
    }else if(method=='mixture'){
        files = unlist(strsplit(file,";"))
        for(i in 1:length(files)){
            fnum=i
            data=read.table(files[i], header=TRUE)
            chrm = data$chromosome[1]
            q25=quantile(data$exp_count, 0.25)
            mf <- model.frame(formula=formula, data=data)
	    mfE <- model.frame(formula=formulaE, data=data)
            X <- model.matrix(attr(mf, "terms"), data=mf)
	    XE <- model.matrix(attr(mfE, "terms"), data=mfE)
	    XNB=as.data.frame(X[,-c(1)])
	    XNBE=as.data.frame(XE[,-c(1)])
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
                    mu1=parms$start$count1
                    mu2=parms$start$count2
                    prob0 =parms$start$zero
                    prop1=parms$prop1
                    prop2=parms$prop2
                    theta1=parms$start$theta1
                    theta2=parms$start$theta2
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
            n <- length(Y)
            kx <- NCOL(X)
            Y0 <- Y <= 0
            Y1 <- Y > 0
            linkstr <- 'logit'
            linkobj <- make.link(linkstr)
            linkinv <- linkobj$linkinv

            #starting params for ze0 component
            model_zero <-.C("pglm_fit", family=as.integer(1), N=as.integer(length(Y)), M=as.integer(ncol(XNB)), y=as.double(Y0), prior=as.double(rep(1,n)), offset=as.double(rep(0,n)), X=as.double(unlist(XNB)),  stratum=as.integer(rep(1,n)),init=as.integer(1), rank=integer(1), Xb=double(n*ncol(XNB)), fitted=as.double((rep(1,n) * Y0 + 0.5)/(rep(1,n) + 1)), resid=double(n), weights=double(n),scale=double(1), df_resid=integer(1), theta=as.double(-1), package='zinba')
            prop0=sum(model_zero$fitted)/n


            #INITIALIZATION OF THE MIXTURE MODEL
           if(initmethod=='quantile'){
	 	if(i == 1){
	                    startprop=startenrichment(c(.15, .001), data, formula, formulaE, initmethod)
	            }
	        data2=data
	        if(sum(colnames(data)=='input_count')==1){data2$input_count=exp(data2$input_count)-1}
		if(sum(colnames(data)=='exp_cnvwin_log')==1){data2$exp_cnvwin_log=exp(data2$exp_cnvwin_log)-1}
		prop2=startprop
            	prop1=1-prop0-prop2
		t=rq(formula, tau=1-prop2, data=data2, method='pfn')
		priorCOUNTweight=rep(10^-10, length(Y))      
		priorCOUNTweight[as.double(which(t$residuals>quantile(t$residuals,1-prop2)))]=1-10^-10
		rm(data2)
  	   }else if(initmethod=='count'){
		    if(i == 1){
	                    prop2=startenrichment(c(.15, .001), data, formula, formulaE,initmethod)
	            }
	            prop1=1-prop0-prop2
	            n1  = round(length(Y) * (1 - prop2))
		    priorCOUNTweight=rep(1-10^-10, length(Y))
                    odY = order(Y)
		    priorCOUNTweight[odY[1:n1]]=10^-10
  	  }else if(initmethod=='pscl'){
            	    mf2 <- model.frame(formula=formula, data=data)
            	    X2 <- model.matrix(attr(mf2, "terms"), data=mf2)
            	    if(i == 1){
                    	a=zeroinfl(formula, data=data,dist='negbin', EM=TRUE)
	            }else{
                        a=zeroinfl(formula, data=data,dist='negbin', EM=TRUE,start=param)
            	    }
            	    leverage=hat(X2, intercept=FALSE)
            	    fdrlevel=threshold
            	    standardized=residuals(a)/sqrt(1-leverage)
            	    pval=1-pnorm(as.matrix(standardized))
            	    fdr=qvalue(pval)
            	    peaks=which(fdr[[3]]<fdrlevel)
		    prop2=length(peaks)/length(Y)
		    prop1=1-prop0-prop2
	            param=list(count=a$coefficients$count, zero=a$coefficients$zero, theta=a$theta)
		    priorCOUNTweight=rep(10^-10, length(Y))      
		    priorCOUNTweight[peaks]=1-10^-10
		    rm(X2)
		    rm(mf2)
		    rm(standardized)
		    rm(leverage)
		    rm(pval)
		    rm(fdr)
		    rm(a)
		    gc()
    	  }
            

	    model_count1 <- .C("pglm_fit", family=as.integer(2), N=as.integer(length(Y)), M=as.integer(ncol(XNB)), y=as.double(Y), prior=as.double(1-priorCOUNTweight), offset=as.double(rep(0,length(Y))), X=as.double(unlist(XNB)),  stratum=as.integer(rep(1,length(Y))),init=as.integer(1), rank=integer(1), Xb=double(length(Y)*ncol(XNB)), fitted=as.double(Y+(Y==0)/6), resid=double(length(Y)), weights=double(length(Y)),scale=double(1), df_resid=integer(1), theta=as.double(-1), package='zinba')
            model_count2 <- .C("pglm_fit", family=as.integer(2), N=as.integer(length(Y)), M=as.integer(ncol(XNBE)), y=as.double(Y), prior=as.double(priorCOUNTweight), offset=as.double(rep(0,length(Y))), X=as.double(unlist(XNB)),  stratum=as.integer(rep(1,length(Y))),init=as.integer(1), rank=integer(1), Xb=double(length(Y)*ncol(XNB)), fitted=as.double(Y+(Y==0)/6), resid=double(length(Y)), weights=double(length(Y)),scale=double(1), df_resid=integer(1), theta=as.double(-1), package='zinba')
            
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
            
            while(abs((ll_old - ll_new)/ll_old) > tol) {
#		    print(ll_new)
                ll_old <- ll[max(1, i-10)]
                prop1=sum(probi1)/n
                prop2=sum(probi2)/n
                #updated values for parameters of component means
                 
                model_zero <- .C("pglm_fit", family=as.integer(1), N=as.integer(length(Y)), M=as.integer(ncol(XNB)), y=as.double(probi0), prior=as.double(rep(1,n)), offset=as.double(rep(0,n)), X=as.double(unlist(XNB)),  stratum=as.integer(rep(1,n)),init=as.integer(1), rank=integer(1), Xb=double(n*ncol(XNB)), fitted=as.double((rep(1,n) * probi0 + 0.5)/(rep(1,n) + 1)), resid=double(n), weights=double(n),scale=double(1), df_resid=integer(1), theta=as.double(-1), package='zinba')  
                model_count1 <- .C("pglm_fit", family=as.integer(0), N=as.integer(length(Y)), M=as.integer(ncol(XNB)), y=as.double(Y), prior=as.double(probi1), offset=as.double(rep(0,length(Y))), X=as.double(unlist(XNB)),  stratum=as.integer(rep(1,length(Y))),init=as.integer(1), rank=integer(1), Xb=double(length(Y)*ncol(XNB)), fitted=as.double(start$count1), resid=double(length(Y)), weights=double(length(Y)),scale=double(1), df_resid=integer(1), theta=as.double(start$theta1), package='zinba')  
                model_count2 <- .C("pglm_fit", family=as.integer(0), N=as.integer(length(Y)), M=as.integer(ncol(XNBE)), y=as.double(Y), prior=as.double(probi2), offset=as.double(rep(0,length(Y))), X=as.double(unlist(XNBE)),  stratum=as.integer(rep(1,length(Y))),init=as.integer(1), rank=integer(1), Xb=double(length(Y)*ncol(XNBE)), fitted=as.double(start$count2), resid=double(length(Y)), weights=double(length(Y)),scale=double(1), df_resid=integer(1), theta=as.double(start$theta2), package='zinba')  

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
		if(i>300){break}
            }
            numpeaks=length(which(probi2>peakconfidence))
	    if(printFullOut == 1){
		data=cbind(data,((data$exp_count>q25)^2),probi2)
		colnames(data)[c(dim(data)[2]-1,dim(data)[2])]=c('q25','peakprob')
	    }else{
	    	data=cbind(data[1:3],((data$exp_count>q25)^2),probi2)
		colnames(data)=c('chromosome','start','stop','q25','peakprob')
	    }
	    if(diff==0){
		data$q25[Y-mui1<0]=0
            }
            line0 = paste("Processing ",files[fnum])
            line2=paste('Selected number of peaks: ', as.character(numpeaks),sep='')
    ### PRINT SIGNIFICANT WINDOWS
            print(paste(c(line0,line2)))
            winfile = paste(winout,"_",chrm,".wins",sep="")
            if(file.exists(winfile)){
                write.table(data,winfile,quote=F,sep="\t",row.names=F,col.names=F,append=T)
            }else{
                write.table(data,winfile,quote=F,sep="\t",row.names=F)
            }
		rm(data); rm(Y); rm(X); rm(XNB); rm(XE);rm(XNBE);rm(probi0); rm(probi1); rm(probi2); rm(mui1); rm(mui2); rm(start); rm(prop1); rm(prop0);gc();
        }
        time.end <- Sys.time()
        print(difftime(time.end,time.start))
        return(winfile)
    }
}
