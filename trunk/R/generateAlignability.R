generateAlignability=function(mapdir, outdir="", athresh=1, extension=0, twoBitFile){
	if(!file.exists(mapdir)){stop("Specified directory 'mapdir' doesnt exist or is not correct")
	}else{
		#if alignability path does not end in /, then put it in
		if( strsplit(mapdir,"")[[1]][length(strsplit(mapdir,"")[[1]])]!='/') mapdir=paste(mapdir, '/', sep='')
	} 
	if(!file.exists(outdir)){
		stop("Specified directory 'outdir' doesnt exist or is not correct")
	}else{
		#if alignability path does not end in /, then put it in
		if( strsplit(outdir,"")[[1]][length(strsplit(outdir,"")[[1]])]!='/') outdir=paste(outdir, '/', sep='')
	} 	 	
	if(!file.exists(twoBitFile)){stop("Specified twoBit file doesnt exist or is not correct")} 
	
	currdir=getwd()
	setwd(mapdir)
	convertmappability(inputfile='map.list', outputfile='temp.wig')
	setwd(currdir)
	alignAdjust(inputfile=paste(mapdir,'temp.wig',sep=''), outDir=outdir,twoBitFile=twoBitFile,athresh=athresh,adjustsize=round(extension/2))
	unlink('temp.wig')
}
