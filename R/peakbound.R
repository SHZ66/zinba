peakbound=function(bpprofile,output,winoffset){

#######################
# USE winoffset for offset size
#######################

    bpProfiles=read.table(bpprofile,header=TRUE,sep="\t")
    bpVector=bpProfiles[,6:ncol(bpProfiles)]
    peakbound2=function(x){
    	##assuming 500 bp window and flanking 500 bp windows
	##TO DO: 
    	##2)make to possibly handle multiple peaks
    	##############################
    	##############################
    	#left bound
    	#find max of center window
    	max=which.max(x[(winoffset+1):(length(x)-winoffset)])+winoffset
    	hold=matrix(0, 1,length(x)-2*winoffset)
    	for(bound in (length(x)-2*winoffset):50){
        	X=(max-bound):max
        	Y=x[(max(max-bound,1)):max]
        	B=sum((X-mean(X))*Y)/sum((X-mean(X))^2)
        	hold[bound]=(B^2)*sum((X-mean(X))^2)/sum((Y-mean(Y))^2)
    	}

    	fit=which.max(hold)+1
    	peakstart=max-fit
    	##############################
    	#right bound
    	hold=matrix(0, 1,length(x)-2*winoffset)
    	for(bound in 50:(length(x)-2*winoffset)){
        	X=max:(max+bound)
        	Y=x[(max:(max+bound))]
        	B=sum((X-mean(X))*Y)/sum((X-mean(X))^2)
        	hold[bound]=(B^2)*sum((X-mean(X))^2)/sum((Y-mean(Y))^2)
    	}
    	fit=which.max(hold)
    	peakend=max+fit
    	##############################
	if(max==(winoffset+1) || max==(length(x)-winoffset)){
		return(c('NA', 'NA','NA', 'NA'))	
	}else{	
	    	return(c(peakstart, peakend, max(x[(winoffset+1):(length(x)-winoffset)]), which.max(x[(winoffset+1):(length(x)-winoffset)])+winoffset))
    	}
    }
    peakCoords=apply(bpVector,1,peakbound2)
    refPeaks1=cbind(bpProfiles[,1:5],t(peakCoords))
    removeNAs=which(refPeaks1[,6]=='NA')
    if(length(removeNAs)>0){
	    refPeaks=refPeaks1[-which(refPeaks1[,6]=='NA'),]
    }else{
	    refPeaks=refPeaks1
    }	
    write.table(refPeaks,output,quote=F,sep="\t",row.names=F,col.names=F, append=TRUE)
}

