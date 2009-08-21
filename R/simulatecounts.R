simulatecounts=function(a, X, numpeaks, peakcount,method){
	if(method=='zicounts'){	
		X=as.matrix(cbind(rep(1, dim(X)[1]),X))
		numcoef=length(a$coefficients)-1
		if(dim(X)[2]!=.5*numcoef){
			warning('number of coefficients do not equal number of columns, check if column for intercept in X')
		}	
		background=exp(X%*%a$coefficients[1:(numcoef/2)])*(1-1/(exp(-X%*%a$coefficients[(numcoef/2+1):numcoef])+1))
		theta=a$coefficients[length(a$coefficients)]
		probz=1/(exp(-X%*%a$coefficients[(numcoef/2+1):numcoef]))
	}
	else if(method=='pscl'){
		X=as.matrix(cbind(rep(1, dim(X)[1]),X))
		numcoef=length(a$coefficients[[1]]) #need to make this more general for diff X and Z
		if(dim(X)[2]!=numcoef){
			warning('number of coefficients do not equal number of columns, check if column for intercept in X')
		}	
		background=exp(X%*%a$coefficients[[1]])*(1-1/(exp(-X%*%a$coefficients[[2]])+1))
		theta=a$theta
		probz=1/(exp(-X%*%a$coefficients[[2]]))	
	}	
	p=numpeaks/length(background)
	if(numpeaks>0){
		enrich=rbinom(length(background),1,p)
	}else{
		enrich=rep(0, length(background))
	}
	count=rep(0, length(background))
	rzinbinom <- function(n, mu, size, zprob) {
		ifelse(runif(n) < zprob, 0, rnbinom(n, mu=mu, size=size))
	}
	for(i in 1:length(background)){
		if(enrich[i]==0){
			count[i]=rzinbinom(1, mu=background[i], size=theta, zprob=probz[i])
			#startpos=matrix((i-1)*500, count, 1)
			#midpoints=round(runif(0,500))
			#write.table(c(midpoints-round(b/2), startpoints+round(b/2)),file, append=TRUE)
		}else{
			count[i]=rzinbinom(1, mu=background[i], size=theta, zprob=probz[i])+rnbinom(1, mu=peakcount, size=theta)
			#startpos=matrix((i-1)*500, count, 1)
			#midpoints=round(runif(0,500))
			#write.table(c(midpoints-round(b/2), startpoints+round(b/2)),file, append=TRUE)
		}
	
	}
	results=as.matrix(cbind(enrich, count))
	colnames(results)=c('enriched', 'count')
	return(results)
}
