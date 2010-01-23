basealigncount=function(inputfile,outputfile,twoBitFile,extension=NULL,filetype="bowtie"){
    if(!file.exists(inputfile)){
    	stop(paste("Input file not found,",inputfile,sep=" "))
    }
    if(!file.exists(twoBitFile)){
        stop(paste("twoBit file not found,",twoBit,sep=" "))
    }
    if(filetype != "bowtie" && filetype != "tagAlign" && filetype != "bed"){
        stop(paste("Incorrect filetype:",filetype,"[bowtie|tagAlign|bed]",sep=" "))
    }
    if(is.null(extension)){
        stop(paste("Need to specify an extension length",,sep=" "))
    }
    cReturn <- .C("baseAlignCounts",as.character(inputfile), as.character(outputfile), as.character(twoBitFile),as.integer(extension),as.character(filetype),PACKAGE="zinba")
}
