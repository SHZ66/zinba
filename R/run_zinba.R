run.zinba=function(seq,align=NULL,input=NULL,cnvarray=NULL,twoBit=NULL,winSize=500,offset=0,aThresh=1,gb=NULL){
    #buildwindowdata=function(seq,align=NULL,input=NULL,cnvarray=NULL,twoBit=NULL,winSize=500,offset=0,aThresh=1,gb=NULL)

#getsigwindows=function(file,covnames,threshold=.01,output)
#basecountimport=function(inputfile, coordfile, outputfile, chromosome='all')

#bpProfiles=read.table(out,header=TRUE,sep="\t")
#bpVector=bpProfiles[,6:(ncol(bpProfiles)-1)]
#peakCoords=apply(bpVector,1,peakbound)
#refPeaks=cbind(bpProfiles[,1:5],t(peakCoords))
#out <- paste(coordFile,"_REFINE_COORDS.txt",sep="")
#write.table(refPeaks,out,quote=F,sep="\t",row.names=F,col.names=F)


#peakbound=function(x)



#time.start <- Sys.time()
#time.end <- Sys.time()
#print(difftime(time.end,time.start))

	Perl.Path <- file.path(.path.package("zimba"), "Perl")
	Fn.Path <- file.path(Perl.Path, "runGetPeaks.pl")
        CMD=paste("--seq", seq,"--align", align,"--gdna", gdna,"--cnvarray", cnvarray,sep=" ")
	system(CMD)
}