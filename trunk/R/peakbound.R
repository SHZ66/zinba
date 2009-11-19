##with new parameters
peakbound=function(bpprofile,output,winSize=200, quantile=.75){
	    localMaximum=function(x, winSize, quantile) {
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
	    	localMax[which(x<quantile(x,quantile, na.rm=TRUE))]=0
	    	#ensure global max is always tagged
	    	localMax[which.max(x)]=1
	    	return(which(localMax>0))
            }


	    #getBPcount <- readLines("MCF7_FAIRE_w500_o125_mix99_InputAllCovs.bpcount",n=-1)
	    #dataLine <- strsplit(x=getBPcount[2],split="\t")

	    maxNumCol = max(count.fields(bpprofile,sep="\t"))
	    bpProfiles=read.table(bpprofile,header=F,col.names=as.character(1:maxNumCol),skip=1,fill=T,sep="\t")
	    searchlength=500
	    peakbound2=function(xx){
		    xna <- is.na(xx) | is.nan(xx) | is.infinite(xx)
		    xc <- xx[!xna]
		    x=as.numeric(xc[6:length(xc)])
		    maxvec=localMaximum(x, winSize=winSize, quantile=quantile)
		    #maxvec=which.max(x)
		    xcoords=matrix(0,length(maxvec),9)
		    bounds=matrix(.C('peakboundc', as.vector(x, mode='integer'), as.integer(length(x)), as.vector(maxvec, mode='integer'), as.integer(length(maxvec)), as.integer(500),  vector("integer",2*length(maxvec)), PACKAGE="zinba")[[6]],length(maxvec),2, byrow=T)
		    xcoords=cbind(matrix(rep(as.matrix(xx[1:5]),length(maxvec)),length(maxvec),5, byrow=T),bounds[,1], bounds[,2],x[maxvec], maxvec)
	   	    xcoords
		    #c(as.matrix(xx[1:5]),bounds[1], bounds[2], x[maxvec], maxvec) 
		}
    refPeaks=do.call('rbind', apply(bpProfiles, 1, peakbound2))
    refPeaks[,c(6,7,9)]=as.numeric(refPeaks[,c(6,7,9)])+as.numeric(refPeaks[,3])-1
    colnames(refPeaks)=c('peakID','chrom','start','stop','strand','peak_start','peak_stop','max_score','max_position')
    refPeaks=as.data.frame(refPeaks)
    if(file.exists(output)){
        write.table(refPeaks,output,quote=F,sep="\t",row.names=F,col.names=F,append=TRUE)
    }else{
        write.table(refPeaks,output,quote=F,sep="\t",row.names=F)
    }
}


