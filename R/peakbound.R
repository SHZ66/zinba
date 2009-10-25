localMaximum=function(x, winSize = 50) 
{
    #from MassSpecWavelet Package in BioConductor, written by Pan Du and Simon Lin, modified by Naim Rashid
    len <- length(x)
    rNum <- ceiling(len/winSize)
    y <- matrix(c(x, rep(x[len], rNum * winSize - len)), nrow = winSize)
    y.maxInd <- apply(y, 2, which.max)
    selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > 
        x[winSize]))
    localMax <- rep(0, len)
    localMax[(selInd - 1) * winSize + y.maxInd[selInd]] <- 1
    shift <- floor(winSize/2)
    rNum <- ceiling((len + shift)/winSize)
    y <- matrix(c(rep(x[1], shift), x, rep(x[len], rNum * winSize - 
        len - shift)), nrow = winSize)
    y.maxInd <- apply(y, 2, which.max)
    selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > 
        x[winSize]))
    localMax[(selInd - 1) * winSize + y.maxInd[selInd] - shift] <- 1
    maxInd <- which(localMax > 0)
    selInd <- which(diff(maxInd) < winSize)
    if (length(selInd) > 0) {
        selMaxInd1 <- maxInd[selInd]
        selMaxInd2 <- maxInd[selInd + 1]
        temp <- x[selMaxInd1] - x[selMaxInd2]
        localMax[selMaxInd1[temp <= 0]] <- 0
        localMax[selMaxInd2[temp > 0]] <- 0
    }
    #delete maxes lower than median
    localMax[which(x<median(x))]=0
    #ensure global max is always tagged
    localMax[which.max(x)]=1
    return(which(localMax>0))
}

peakbound=function(bpprofile,output,winoffset=0){
	    bpProfiles=read.table(bpprofile,header=TRUE,sep="\t")
	    searchlength=500
	    peakbound2=function(xx){	
		    x=as.numeric(xx[6:length(xx)])
		    #maxvec=localMaximum(x)
		    maxvec=which.max(x)
		    xcoords=matrix(0,length(maxvec),9)
		    start=as.numeric(xx[3])-1
         	    bounds=matrix(.C('peakboundc', as.vector(x, mode='integer'), as.integer(length(x)), as.vector(maxvec, mode='integer'), as.integer(length(maxvec)), as.integer(500),  vector("integer",2*length(maxvec)), PACKAGE="basecount")[[6]],length(maxvec),2, byrow=T)+start
			xcoords=cbind(matrix(rep(as.matrix(xx[1:5]),length(maxvec)),length(maxvec),5, byrow=T),bounds, x[maxvec], maxvec+start)
		    		
	   colnames(xcoords)=c(names(xx[1:5]), c('peak_start','peak_stop','max_score','max_position'))
	   xcoords
		}

    	     		
	   
    refPeaks=as.data.frame(do.call('rbind', apply(bpProfiles, 1, peakbound2)))
    refPeaks[,c(6,7,9)]=refPeaks[,c(6,7,9)]+refPeaks[,3]-1
    if(file.exists(output)){
        write.table(refPeaks,output,quote=F,sep="\t",row.names=F,col.names=F,append=TRUE)
    }else{
        write.table(refPeaks,output,quote=F,sep="\t",row.names=F)
    }
}

