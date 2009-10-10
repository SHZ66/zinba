peakbound=function(bpprofile,output,winoffset=0){
    bpProfiles=read.table(bpprofile,header=TRUE,sep="\t")
    bpVector=as.matrix(bpProfiles[,6:ncol(bpProfiles)])
    peakbound2=function(x){
	##assuming 500 bp window and flanking 500 bp windows
	##TO DO: 
	##2)make to possibly handle multiple peaks
	##############################
	##############################
	#left bound
	#find max of center window
	searchlength=(length(x)-1)/2
	max=searchlength+1
	hold=matrix(0, 1,searchlength)
	for(bound in searchlength:50){
	    X=(max-bound):max
	    Y=x[(max(max-bound,1)):max]
	    B=sum((X-mean(X))*Y)/sum((X-mean(X))^2)
	    hold[bound]=(B^2)*sum((X-mean(X))^2)/sum((Y-mean(Y))^2)
	}
	fit=which.max(hold)+1
	peakstart=max-fit
	##############################
	#right bound
	hold=matrix(0, 1,searchlength)
	for(bound in 50:searchlength){
	    X=max:(max+bound)
	    Y=x[max:(max(max+bound,searchlength))]
	    B=sum((X-mean(X))*Y)/sum((X-mean(X))^2)
	    hold[bound]=(B^2)*sum((X-mean(X))^2)/sum((Y-mean(Y))^2)
	}
	fit=which.max(hold)
	peakend=max+fit
	##############################

        return(c(peakstart, peakend,max(x),which.max(x)))
    }
    peakCoords=apply(bpVector,1,peakbound2)
    refPeaks=cbind(bpProfiles[,1:5],t(peakCoords))
    refPeaks[,c(6,7,9)]=refPeaks[,c(6,7,9)]+refPeaks[,3]-1

#   refPeaks1=cbind(bpProfiles[,1:5],t(peakCoords))
#   removeNAs=which(refPeaks1[,6]=='NA')
#   if(length(removeNAs)>0){
#	refPeaks=refPeaks1[-which(refPeaks1[,6]=='NA'),]
#	refPeaks[,c(6,7,9)]=refPeaks[,c(6,7,9)]+refPeaks[,3]-1
#    }else{
#	refPeaks=refPeaks1
#	refPeaks[,c(6,7,9)]=refPeaks[,c(6,7,9)]+refPeaks[,3]-1
#    }
#    refPeaks=refPeaks[order(refPeaks[,9]),]
#    #number of unique maxes
#    uniquemaxes=length(levels(factor(refPeaks[,9])))
#    maxes=as.numeric(levels(factor(refPeaks[,9])))
#    for(i in 1:uniquemaxes){
#	#gets rows corresponding to unique max i		
#	maxmatch=which(as.numeric(refPeaks[,9])==maxes[i])
#	if(length(maxmatch)>1){
#	    minstart=min(refPeaks[maxmatch,6]) #min of all start for those with same max
#	    maxend=max(refPeaks[maxmatch,7]) #max of all end for those with same max
#	    #save information for first corresponding row to new matrix
#	    refPeaks[maxmatch,6]=minstart
#	    refPeaks[maxmatch,7]=maxend
#	    #delete duplicates
#	    refPeaks=refPeaks[-maxmatch[2:length(maxmatch)],]
#	    #replace start and stop bounds of row with new merged start and stop
#	}
#   }

    if(file.exists(output)){
        write.table(refPeaks,output,quote=F,sep="\t",row.names=F,col.names=F,append=TRUE)
    }else{
        write.table(refPeaks,output,quote=F,sep="\t",row.names=F)
    }
#    write.table(refPeaks,output,quote=F,sep="\t",row.names=F,col.names=F,append=TRUE)
}





