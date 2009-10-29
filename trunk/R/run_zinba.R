#run.zinba=function(paramFile=NULL,seq=NULL,align=NULL,input=NULL,cnvarray=NULL,cnvexpwin=NULL,cnvexpcustom=NULL,rand=NULL,twoBit=NULL,winSize=500,offset=0,percN=0.1,inputTrans=NULL,inputRandLog=NULL,gb="hg18",basecountfile=NULL,threshold=0.01,numProc=1,buildwin=0,winoffset=0,refinepeaks=1,printLog=0,method='pscl'){
run.zinba=function(paramFile=NULL,seq=NULL,align=NULL,input=NULL,twoBit=NULL,winSize=500,offset=0,cnvWinSize=100000,cnvOffset=0,basecountfile=NULL,threshold=0.01,numProc=1,buildwin=0,winoffset=500,refinepeaks=1,printLog=0,method='pscl'){
	time.start <- Sys.time()
        if(buildwin==1){
            buildwindowdata(paramFile=paramFile,seq=seq,align=align,input=input,twoBit=twoBit,winSize=winSize,offset=offset,cnvWinSize=cnvWinSize,cnvOffset=cnvOffset)
	}
	Perl.Path <- file.path(.path.package("zinba"), "exec")
	Fn.Path <- file.path(Perl.Path, "runGetPeaks_v2.pl")
	if(refinepeaks==1 && is.null(basecountfile)){
		print(paste("Basecount file must be specified, currently",basecountfile,sep=" "))
	}else if (is.null(paramFile)){
		print(paste("Need parameter file with formula specified ",paramFile,sep=" "))
	}else{
                print(paste("Getting significant windows"))
                print(paste("Using parameter file",paramFile,sep=" "))
                print(paste("Running",numProc,"parallel jobs",sep=" "))
	        CMD=paste(Fn.Path,"--param-file", paramFile,"--threshold",threshold,"--win-size",winSize,"--basecount_file",basecountfile,"--method",method,"--processes",numProc,"--refine_peaks",refinepeaks,"--win-offset",winoffset,sep=" ")
		if(printLog==1){
			CMD=paste(CMD,"--print-log 1",sep=" ")
		}
		system(CMD)
	}
	time.end <- Sys.time()
	print(difftime(time.end,time.start))
}