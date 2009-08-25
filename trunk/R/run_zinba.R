run.zinba=function(win_list=NULL,seq,align=NULL,input=NULL,cnvarray=NULL,twoBit=NULL,winSize=500,offset=0,aThresh=1,percN=0.1,gb=NULL,basecountfile=NULL,covs=NULL,threshold=0.01,numProc=1,method='zicounts'){
	time.start <- Sys.time()
        if(is.null(win_list)){
            buildwindowdata(seq=seq,align=align,input=input,cnvarray=cnvarray,twoBit=twoBit,winSize=winSize,offset=offset,aThresh=aThresh,percN=percN,gb=gb,nProcess=numProc)
            data_list <- read.table(paste(seq,".list",sep=""))
	}else{
            data_list <- read.table(list)
        }
        

#######
	Perl.Path <- file.path(.path.package("zinba"), "exec")
	Fn.Path <- file.path(Perl.Path, "runGetPeaks_v2.pl")
        CMD=paste(Fn.Path,"--win-file", data_list,"--threshold",threshold,"--win-size",winSize,"--basecount_file",basecountfile,"--method",method,"--covs",covs,"--processes",numProc,sep=" ")
	system(CMD)
######
        
#        for (i in 1:nrow(data_list)){
#	    print(paste("Getting windows for ",as.character(data_list[i,1])))
#            coordout <- paste(substr(data_list[i,1],1,((nchar(as.character(data_list[i,1])))-4)),"_PEAK_COORDS.temp",sep="")
#            winout <- paste(substr(data_list[i,1],1,((nchar(as.character(data_list[i,1])))-4)),".wins",sep="")
#            peakout <- paste(substr(data_list[i,1],1,((nchar(as.character(data_list[i,1])))-4)),".peaks",sep="")
#            bpOut <- paste(coordout,"_BPcount",sep="")
#
#	    getsigwindows(file=as.character(data_list[i,1]),covnames=covs,threshold=threshold,winout=winout,coordout=coordout,offset=(winSize/2),method=method)
#	    basecountimport(inputfile=basecountfile,coordfile=coordout,outputfile=bpOut,chromosome=data_list[i,2])
#	    peakbound(profile=bpOut,output=peakout)
#
#            unlink(coordout)
#            unlink(bpOut)
#	    gc()
#        }
        unlink(paste(seq,".list",sep=""))
	time.end <- Sys.time()
	print(difftime(time.end,time.start))
}