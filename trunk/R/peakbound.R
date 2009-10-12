peakbound=function(bpprofile,output,winoffset=0){
    bpProfiles=read.table(bpprofile,header=TRUE,sep="\t")
    bpVector=as.matrix(bpProfiles[,6:ncol(bpProfiles)])
    x=bpVector
	searchlength=(dim(x)[2]-1)/2
	max=searchlength+1
	hold=matrix(0, dim(x)[1],searchlength)
	for(bound in searchlength:50){
	    X=(max-bound):max
	    meanx=mean(X)
	    Y=x[,(max(max-bound,1)):max]
	    meany=apply(Y, 1, mean)
	    B=(Y%*%(X-mean(X)))/sum((X-mean(X))^2)
	    hold[,bound]=(B^2)*sum((X-meanx)^2)/(((Y-meany)^2)%*%rep(1,dim(Y)[2]))
	}
	fit=apply(hold, 1, which.max)
	peakstart=max-fit
	##############################
	#right bound
	hold=matrix(0, dim(x)[1],searchlength)
	for(bound in 50:searchlength){
	    X=max:(max+bound)
	    meanx=mean(X)
	    Y=x[,max:(max(max+bound,searchlength))]
   	    meany=apply(Y, 1, mean)
	    B=(Y%*%(X-mean(X)))/sum((X-mean(X))^2)
	    hold[,bound]=(B^2)*sum((X-meanx)^2)/(((Y-meany)^2)%*%rep(1,dim(Y)[2]))
	}
	fit=apply(hold, 1, which.max)
	peakend=max+fit
	##############################
    refPeaks=cbind(bpProfiles[,1:5],peakstart, peakend,x[,251],rep(251, dim(x)[1]))
    refPeaks[,c(6,7,9)]=refPeaks[,c(6,7,9)]+refPeaks[,3]-1

    colnames(refPeaks)[c(dim(refPeaks)[2]-3,dim(refPeaks)[2]-2,dim(refPeaks)[2]-1, dim(refPeaks)[2])]=c('peak_start','peak_stop','max_score','max_position')
    if(file.exists(output)){
        write.table(refPeaks,output,quote=F,sep="\t",row.names=F,col.names=F,append=TRUE)
    }else{
        write.table(refPeaks,output,quote=F,sep="\t",row.names=F)
    }
}

