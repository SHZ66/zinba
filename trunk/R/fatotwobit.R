fatotwobit=function(faFile=" ",fadir = " ",outFile=NULL){
    	if(fadir!=" " & !file.exists(fadir)){
		stop("Specified directory 'fadir' doesnt exist or is not correct")
	}else{
		if(strsplit(fadir,"")[[1]][length(strsplit(fadir,"")[[1]])]!='/' & fadir!=" ") fadir=paste(fadir, '/', sep='')
	}
   
	if(fadir!=" " & faFile==" "){
		fafiles=dir(fadir, pattern="\\.fa")
		for(i in 1:length(fafiles)){
			print(paste("processing",paste(fadir,fafiles[i],sep="")))
			.C("faToTwoBit",as.character(paste(fadir,fafiles[i],sep="")),as.character(outFile),PACKAGE="zinba")
		}
	}else if(faFile!=" " & fadir==" "){
		.C("faToTwoBit",as.character(faFile),as.character(outFile),PACKAGE="zinba")
	}else{
		stop("Either need to specify a directory containing .fa files to be converted or specify a specific faFile")
	}
	
}

