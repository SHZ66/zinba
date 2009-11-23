run.zinba=function(filelist=NULL,formula=NULL,formulaE=NULL,outfile=NULL,seq=NULL,align=NULL,input="none",twoBit=NULL,winSize=500,offset=0,cnvWinSize=100000,cnvOffset=0,basecountfile=NULL,threshold=0.01,peakconfidence=.8,tol=10^-5,numProc=1,buildwin=0,pWinSize=200,pquant=0.75,refinepeaks=1,printLog=0,method='pscl',initmethod='count'){
        library(multicore)
        library(doMC)
        library(foreach)
	time.start <- Sys.time()
	if(is.null(formulaE)){
		formulaE=formula
	}
        if(buildwin==1){
            buildwindowdata(seq=seq,align=align,input=input,twoBit=twoBit,winSize=winSize,offset=offset,cnvWinSize=cnvWinSize,cnvOffset=cnvOffset,filelist=filelist)
	}
	if(refinepeaks==1 && is.null(basecountfile)){
		stop(paste("Basecount file must be specified, currently",basecountfile,sep=" "))
	}else if (is.null(filelist)){
		stop(paste("Need list of files ",filelist,sep=" "))
	}else if(method != 'pscl' && method != 'mixture'){
		stop(paste("Method should be either pscl or mixture, currently",method))
	}else{
            params=scan(filelist,what=character(0))
#            winout=paste(outfile,".wins",sep="")
	    winlist=paste(outfile,".winlist",sep="")
            peakout=paste(outfile,".peaks",sep="")
            coordout=paste(outfile,".coords",sep="")
            bpout=paste(outfile,".bpcount",sep="")

            registerDoMC(numProc)
            mcoptions <- list(preschedule = FALSE, set.seed = FALSE)
            getDoParWorkers()
            winfiles <- foreach(i=1:length(params),.combine='rbind',.inorder=FALSE,.errorhandling="remove",.options.multicore = mcoptions) %dopar%
                getsigwindows(file=params[i],formula=formula,threshold=threshold,winout=outfile,peakconfidence=peakconfidence,tol=tol,method=method,initmethod=initmethod )
 
	    write.table(winfiles,winlist,quote=F,row.names=F,col.names=F)
	    #collapsewins(winlist=winlist,winout=winout)
	    if(refinepeaks==1){
		getrefinedpeaks(winlist=winlist,coordout=coordout,basecountfile=basecountfile,bpout=bpout,peakout=peakout,twoBit=twoBit,winSize=winSize,pWinSize=pWinSize,pquant=pquant,peakconfidence=peakconfidence,threshold=threshold,method=method)
	    }
	}
	time.end <- Sys.time()
	print(difftime(time.end,time.start))
}
