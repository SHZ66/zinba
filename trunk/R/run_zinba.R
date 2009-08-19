run.zinba=function(seq,align=NULL,input=NULL,cnvarray=NULL,twoBit=NULL,winSize=500,offset=0,aThresh=1,gb=NULL,basecountfile=NULL){

	time.start <- Sys.time()
	buildwindowdata(seq,align,input,cnvarray,twoBit,winSize,offset,aThresh,gb)
    
    #need to run getsigwins for each chrm and offset and assign output file;
    #this prints 2 files output (.wins) and output_PEAK_COORDS.temp
#	getsigwindows(file,covnames,threshold=.01,output)

    #coordfile is output from previous
#	basecountimport(inputfile, coordfile, outputfile, chromosome='all')

#	peakbound=function(profile,output)

	time.end <- Sys.time()
	print(difftime(time.end,time.start))
}