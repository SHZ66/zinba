run.zinba=function(win_list=NULL,seq=NULL,align=NULL,input=NULL,cnvarray=NULL,cnvexpwin=NULL,cnvexpcustom=NULL,rand=NULL,twoBit=NULL,winSize=500,offset=0,percN=0.1,inputTrans=NULL,inputRandLog=NULL,gb="hg18",basecountfile=NULL,covs="",threshold=0.01,numProc=1,printLog=0,method='pscl'){
	time.start <- Sys.time()
        if(is.null(win_list)){
            buildwindowdata(seq=seq,align=align,input=input,cnvarray=cnvarray,cnvexpwin=cnvexpwin,cnvexpcustom=cnvexpcustom,rand=rand,twoBit=twoBit,winSize=winSize,offset=offset,inputTrans=inputTrans,inputRandLog=inputRandLog,percN=percN,gb=gb,nProcess=numProc)
            win_list <- paste(seq,".list",sep="")
	}
        print(paste("Processing files listed in",win_list,sep=" "))
	Perl.Path <- file.path(.path.package("zinba"), "exec")
	Fn.Path <- file.path(Perl.Path, "runGetPeaks_v2.pl")
	if(is.null(basecountfile)){
		print(paste("Basecount file must be specified, currently",basecountfile,sep=" "))
	}else{
	        CMD=paste(Fn.Path,"--win-file", win_list,"--threshold",threshold,"--win-size",winSize,"--basecount_file",basecountfile,"--method",method,"--covs",covs,"--processes",numProc,sep=" ")
		if(printLog==1){
			CMD=paste(CMD,"--print-log 1",sep=" ")
		}
		system(CMD)
	}
	time.end <- Sys.time()
	print(difftime(time.end,time.start))
}