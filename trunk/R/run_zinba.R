run.zinba=function(win_list=NULL,seq,align=NULL,input=NULL,cnvarray=NULL,twoBit=NULL,winSize=500,offset=0,aThresh=1,percN=0.1,gb=NULL,basecountfile=NULL,covs=NULL,threshold=0.01,numProc=1,method='zicounts'){
	time.start <- Sys.time()
        if(is.null(win_list)){
            buildwindowdata(seq=seq,align=align,input=input,cnvarray=cnvarray,twoBit=twoBit,winSize=winSize,offset=offset,aThresh=aThresh,percN=percN,gb=gb,nProcess=numProc)
            win_list <- paste(seq,".list",sep="")
	}

        print(paste("Processing files listed in",win_list,sep=" "))
	Perl.Path <- file.path(.path.package("zinba"), "exec")
	Fn.Path <- file.path(Perl.Path, "runGetPeaks_v2.pl")
        CMD=paste(Fn.Path,"--win-file", win_list,"--threshold",threshold,"--win-size",winSize,"--basecount_file",basecountfile,"--method",method,"--covs",covs,"--processes",numProc,sep=" ")
	system(CMD)

	time.end <- Sys.time()
	print(difftime(time.end,time.start))
}