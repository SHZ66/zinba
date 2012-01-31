basecountimport=function(inputfile,winlist,threshold=.01,method='pscl',printFullOut=0,outputfile,twobitfile,chromosome='all', winGap=0, FDR=FALSE){
    if(!file.exists(inputfile)){
    	stop(paste("Input file not found,",inputfile,sep=" "))
    }
    if(file.info(inputfile)$size==0){
    	stop(paste("Input file ",inputfile, "has size 0, check for disk space issues",sep=" "))
    }
    if(!file.exists(twobitfile)){
        stop(paste("twoBit file not found,",twobitfile,sep=" "))
    }
    if(!file.exists(winlist)){
        stop(paste("Window list file not found,",winlist,sep=" "))
    }
    if(is.null(outputfile)){
    	stop(paste("Need to specify an outputfile,",outputfile,sep=" "))
    }

		if(length(grep(".bin.gz", inputfile))>0){
			cat("Using compressed binary version of basecount file\n")
			cat("Uncompressing file", inputfile,"\n")
			library(R.utils)
			gunzip(inputfile, remove=F, overwrite=T)
			binary=1
		}else{
			binary=0
		}

    cReturn <- .C("getSeqCountProfile",as.character(inputfile),as.character(winlist),
				as.double(threshold),as.character(method),as.integer(printFullOut),
				as.character(outputfile),as.character(twobitfile),as.character(chromosome), 
				as.integer(winGap), as.integer(FDR^2),as.integer(binary),PACKAGE="zinba")
	
		if(binary==1) unlink(gsub("\\.gz$", "", inputfile))
}
