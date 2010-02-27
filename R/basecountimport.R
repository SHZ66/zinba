basecountimport=function(inputfile,winlist,threshold=.01,method='pscl',printFullOut=0,outputfile,twobitfile,chromosome='all'){
    if(!file.exists(inputfile)){
    	stop(paste("Input file not found,",inputfile,sep=" "))
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
    cReturn <- .C("getSeqCountProfile",as.character(inputfile),as.character(winlist),as.double(threshold),as.character(method),as.integer(printFullOut),as.character(outputfile),as.character(twobitfile),as.character(chromosome),PACKAGE="zinba")
}
