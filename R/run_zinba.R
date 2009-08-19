run.zinba=function(seq,align=NULL,input=NULL,cnvarray=NULL,twoBit=NULL,winSize=500,offset=0,aThresh=1,gb=NULL,basecountfile=NULL,covs=NULL,threshold=0.01){
	time.start <- Sys.time()
	buildwindowdata(seq=seq,align=align,cnvarray=cnvarray,twoBit=twoBit,offset=offset,aThresh=aThresh,gb=gb)    

	data_list <- read.table(paste(seq,".list",sep=""))
        for (i in 1:nrow(data_list)){
            coordout <- paste(substr(data_list[i,1],1,((nchar(as.character(data_list[i,1])))-4)),"_PEAK_COORDS.temp",sep="")
            winout <- paste(substr(data_list[i,1],1,((nchar(as.character(data_list[i,1])))-4)),".wins",sep="")
            peakout <- paste(substr(data_list[i,1],1,((nchar(as.character(data_list[i,1])))-4)),".peaks",sep="")
            bpOut <- paste(coordout,"_BPcount",sep="")

	    getsigwindows(file=data_list[i,1],covnames=covs,threshold=threshold,winout=winout,coordout=coordout,offset=(winSize/2))
	    basecountimport(inputfile=basecountfile,coordfile=coordout,outputfile=bpOut,chromosome=data_list[i,2])
	    peakbound=function(profile=bpOut,output=peakout)

#            unlink(coordout)
#            unlink(winout)
#            unlink(bpOut)
        }
#        unlink(paste(seq,".list",sep=""))
	time.end <- Sys.time()
	print(difftime(time.end,time.start))
}