basealigncount=function(inputfile,outputfile,twoBitFile,extension=NULL,filetype="bowtie", binary=0){
    if(!file.exists(inputfile)){
    	stop(paste("Input file not found,",inputfile,sep=" "))
    }
    if(!file.exists(twoBitFile)){
        stop(paste("twoBit file not found,",twoBitFile,sep=" "))
    }
    if(filetype != "bowtie" && filetype != "tagAlign" && filetype != "bed"){
        stop(paste("Incorrect filetype:",filetype,"[bowtie|tagAlign|bed]",sep=" "))
    }
    if(is.null(extension)){
        stop(paste("Need to specify an extension length",,sep=" "))
    }
    cReturn <- .C("baseAlignCounts",as.character(inputfile), as.character(outputfile), as.character(twoBitFile),as.integer(extension),as.character(filetype),as.integer(binary), PACKAGE="zinba")
		if(binary==1){
			library(R.utils)
			cat("Because binary option is selected, appending '.bin.gz' to output path.\n")
			newpath=paste(outputfile, ".bin.gz", sep="")
			gzip(outputfile, destname=newpath, overwrite=T, remove=TRUE)
			cat("Compressed binary basecount file can be found at", newpath, "\n")
		}
}
