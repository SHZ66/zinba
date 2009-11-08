run.zinba=function(paramFile=NULL,formula=NULL,outfile=NULL,seq=NULL,align=NULL,input=NULL,twoBit=NULL,winSize=500,offset=0,cnvWinSize=100000,cnvOffset=0,basecountfile=NULL,threshold=0.01,peakconfidence=.8,tol=10^-5,numProc=1,buildwin=0,pWinSize=200,pquant=0.75,refinepeaks=1,printLog=0,method='pscl'){
        library(multicore)
        library(doMC)
        library(foreach)
	time.start <- Sys.time()
        if(buildwin==1){
            buildwindowdata(seq=seq,align=align,input=input,twoBit=twoBit,winSize=winSize,offset=offset,cnvWinSize=cnvWinSize,cnvOffset=cnvOffset)
	}
	if(refinepeaks==1 && is.null(basecountfile)){
		print(paste("Basecount file must be specified, currently",basecountfile,sep=" "))
	}else if (is.null(paramFile)){
		print(paste("Need parameter file with formula specified ",paramFile,sep=" "))
	}else{
            params=scan(paramFile,what=character(0))
            winout=paste(outfile,".wins",sep="")
            peakout=paste(outfile,".peaks",sep="")
            coordout=paste(outfile,".coords",sep="")
            bpout=paste(outfile,".bpcount",sep="")

            registerDoMC(numProc)
            mcoptions <- list(preschedule = FALSE, set.seed = FALSE)
            getDoParWorkers()
            rf <- foreach(i=1:length(params),.options.multicore = mcoptions) %dopar%
                getsigwindows(file=params[i],formula=formula,threshold=threshold,winout=winout,peakconfidence=peakconfidence,tol=tol,method=method)

	    if(refinepeaks==1){
		getrefinedpeaks(winout=winout,coordout=coordout,basecountfile=basecountfile,bpout=bpout,peakout=peakout,twoBit=twoBit,pWinSize=pWinSize,pquant=pquant,peakconfidence=peakconfidence,threshold=threshold,method=method)
	    }
	}
	time.end <- Sys.time()
	print(difftime(time.end,time.start))
}