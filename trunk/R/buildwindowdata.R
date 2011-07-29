buildwindowdata=function(seq,input="none",align,twoBit,winSize=500,offset=0,cnvWinSize=100000,cnvOffset=0,filelist,filetype="bowtie",extension, outdir="default"){
	if(!file.exists(align)){
		stop(paste("Specified Alignability directory not found,",align,sep=" "))
	}else{
		#if alignability path does not end in /, then put it in
		if( strsplit(align,"")[[1]][length(strsplit(align,"")[[1]])]!='/') align=paste(align, '/', sep='')
	}	
	if(!file.exists(seq)){
		stop(paste("Seq file not found,",seq,sep=" "))
	}
	if(input != "none" && !file.exists(input)){
		stop(paste("Input file not found,",input,sep=" "))
	}
	if(!file.exists(twoBit)){
		stop(paste("twoBit file not found,",twoBit,sep=" "))
	}
	if(filetype != "bowtie" && filetype != "tagAlign" && filetype != "bed"){
		stop(paste("Incorrect filetype:",filetype,"[bowtie|tagAlign|bed]",sep=" "))
	}
	if(is.null(extension)){
		stop(paste("extension was not specified"))
	}
	if(is.null(filelist)){
		stop(paste("filelist was not specified"))
	}
	if(!file.exists(outdir) && outdir!="default"){
		stop(paste("buildwindowdata: specified output directory to place built datafiles doesnt exist"))
	}
	wigfiles <- file.access(dir(align,pattern="wig",full.names=T))
	twigfiles <- table(wigfiles)
	for(i in 1:length(twigfiles)){
		if(twigfiles[[i]] == -1){
			print(paste("Unable to access some wig files in align directory:",align,sep=" "))
		}
	}
	cReturn <- .C("buildWindows",as.character(seq),as.character(input),as.character(align),as.character(twoBit),as.integer(winSize),as.integer(offset),as.integer(cnvWinSize),as.integer(cnvOffset),as.character(filetype),as.character(filelist),as.integer(extension), as.character(outdir),PACKAGE="zinba")
	gc()
}
