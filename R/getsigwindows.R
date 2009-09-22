#for covnames, enter desired covariates as c('gcPerc', 'align','cnv') for example
getsigwindows=function(file,covnames,threshold=.01,winout,coordout,offset=0, method='zicounts'){
	time.start <- Sys.time()
	library(zicounts)
	library(qvalue)
        library(pscl)
        options(scipen=999)
	###### USER INPUT############################
	data=read.table(file, header=TRUE)
	if(method=='pscl'){
		cov_text=paste(covnames, sep='', collapse='+')
		covariates=eval(parse(text=paste("model.matrix(~", paste(covnames, sep='', collapse='+'), ",data)")))
		modelcommand=paste("zeroinfl(exp_count ~ ",cov_text,",data=data,dist = 'negbin', EM=TRUE)")
	}else if(method=='zicounts'){
		covariates=eval(parse(text=paste("model.matrix(~", paste(covnames, sep='', collapse='+'), ",data)")))
		cov_text=paste(covnames, sep='', collapse='+')
		modelcommand=paste("zicounts(resp = exp_count ~ ., x =  ~ ",cov_text,", z =  ~", cov_text,",data=data)")
	}
	fdrlevel=threshold
	##############################################
	#########PROCESSING STEPS########################
	#################################################
	#LEVERAGE CALCULATION############################
	leverage=hat(covariates, intercept=FALSE)
	########output is 'leverage'#####################
	#################################################
	#RESIDUAL STANDARDIZATION########################
	#TO DO:  unlink design matrices for count and zero models (Z neq X)
	a=eval(parse(text=modelcommand))
        
	if(method=="zicounts"){
        	if(as.character(a$se[length(a$se)])=="NaN" || a$coefficients[length(a$coefficients)]>15){
			#checks for NaN in theta or hugely inflated then, switches to pscl if using zicounts
			print("Switching to zeroinfl")
                  	print(paste("For ",file," se is ",as.character(a$se[length(a$se)]),sep=""))
                        print(paste("For ",file," coeff is ",a$coefficients[length(a$coefficients)],sep=""))			
			cov_text=paste(rep('data$', length(covnames)), covnames, sep='', collapse='+')
			modelcommand=paste("zeroinfl(data$exp_count ~ ",cov_text,",dist = 'negbin', EM=TRUE)")
			a=eval(parse(text=modelcommand))
			#now switch method to pscl for residual estimation
			method='pscl'
                        print(paste("se now ",as.character(a$se[length(a$se)]),sep=""))
                        print(paste("coeff now ",a$coefficients[length(a$coefficients)],sep=""))
		}
	}
#	print(summary(a))
      	if(method=='zicounts'){
		link=make.link('logit')
		linkinv=link$linkinv
		X=covariates
		coefc=a$coefficients[1:dim(X)[2]]
		coefz=a$coefficients[(dim(X)[2]+1):(2*dim(X)[2])]
        	mu <- exp(as.matrix(X) %*% coefc)[,1]
        	phi <- linkinv(as.matrix(X) %*% coefz)[,1]
        	Yhat <- (1-phi) * mu
        	res <- a$data - Yhat
	        theta1 <- 1/exp(a$coefficients[2*(length(covnames)+1)+1])
		vv <- Yhat * (1 + (phi + theta1) * mu)
		standardized=res/sqrt(vv*(1-leverage))
	}else if(method=='pscl'){
		standardized=residuals(a)/sqrt(1-leverage)
	}
	########output is 'standardized'#################
	#################################################
	#THRESHOLDING/PEAK CALLING#######################
	#TO DO:  implement weighting to remove sig small count windows
	pval=1-pnorm(as.matrix(standardized))
	fdr=qvalue(pval)
	numpeaks=length(which(fdr[[3]]<fdrlevel))
	minresid=min(standardized[which(fdr[[3]]<fdrlevel)])
	line1='|Selection Summary|'
	line2=paste('Selected number of peaks: ', as.character(numpeaks), sep='')
	line3=paste('Minimum Standardized Residual Value of peaks: ', as.character(minresid), sep='')
	sigpeaks=cbind(data[which(fdr[[3]]<fdrlevel),], fdr[[3]][which(fdr[[3]]<fdrlevel)], standardized[which(fdr[[3]]<fdrlevel)])
	colnames(sigpeaks)[c(dim(sigpeaks)[2]-1, dim(sigpeaks)[2])]=c('q-value', 'residual')

### FORMAT PEAK COORDINATE DATA
	peakID=paste(sigpeaks$chromosome,sigpeaks$start,sigpeaks$end,sep=":")
	coordinates=cbind(peakID,as.character(sigpeaks$chromosome),(sigpeaks$start-offset),(sigpeaks$end+offset),"+")
	write.table(coordinates,coordout,quote=F,sep="\t",row.names=F,col.names=F)
###############################

#	sigpeaksSORTED=sigpeaks[order(sigpeaks[,dim(data)[2]+2], decreasing=TRUE),]
#	colnames(sigpeaksSORTED)[c(dim(sigpeaksSORTED)[2]-1, dim(sigpeaksSORTED)[2])]=c('q-value', 'residual')
	#######output is 'sigpeaks' and 'sigpeaksSORTED'#
	#################################################
	print(c(line1, line2,line3))
	write.table(sigpeaks,winout,quote=F,sep="\t",row.names=F)
      	time.end <- Sys.time()
	print(difftime(time.end,time.start))
}
