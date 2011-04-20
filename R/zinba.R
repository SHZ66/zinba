zinba=function(outfile=NULL,seq=NULL,align=NULL,input="none",twoBit=NULL,basecountfile=NULL,threshold=0.05,numProc=1,
refinepeaks=1,printFullOut=0,filetype="bowtie",extension=NULL,broad=FALSE){
	if(is.null(outfile)){stop("output prefix must be specified")}
	if(is.null(seq)){stop("path to mapped experimental reads must be specified")}
	if(!file.exists(seq)){stop("sequencing file does not exist at the specified path, check path for errors")}
	if(!file.exists(align)){stop("mappability directory does not exist at the specified path, check path for errors")}
	if(!file.exists(input) & input!="none"){stop("sequencing file does not exist at the specified path, check path for errors")}
	if(!file.exists(twoBit)){stop(".2bit file does not exist at the specified path, check path for errors")}
	if(is.null(basecountfile) & refinepeaks==1){stop("basecount file needs to be specified for peak refinement")}
	if(!file.exists(basecountfile) & refinepeaks==1){stop("basecount file does not exist at the specified path, check path for errors")}
	if(is.null(extension)){stop("sequencing read extension length (average length of fragment library) needs to be specified")}

	buildwin=1
	winSize=250	
	offset=125
	cnvWinSize=100000
	cnvOffset=2500
	pquant=1
	initmethod="count"
	printFullOut=1
	diff=0
	pWinSize=200
	method="mixture"
	selectmodel=TRUE
	selectchr="chr22"
	selecttype="dirty"
	winGap=0
	
	if(broad==TRUE) winGap=5000
	if(input=="none" & broad==FALSE){ 
		selectcovs=c("gcPerc", "align_perc", "exp_cnvwin_log")
	}else if(input=="none" & broad==TRUE){
		selectcovs=c("gcPerc", "align_perc")
	}else if(input!="none"){
		selectcovs=c("input_count")
		selecttype="complete"
	}


	run.zinba(
		align=align,
		numProc=numProc,
		seq=seq,
		basecountfile=basecountfile,
		filetype=filetype,
		offset=125,
		buildwin=1,
		outfile=outfile,
		threshold=threshold,
		twoBit=twoBit,
		cnvOffset=2500,
		pquant=1,
		winGap=winGap,
		cnvWinSize=100000,	
		initmethod="count",
		printFullOut=1,
		winSize=250,
		diff=0,
		pWinSize=200,	
		extension=extension,
		method="mixture",
		refinepeaks=refinepeaks,
		selectmodel=TRUE,
		selectchr="chr22",
		selecttype=selecttype,
		selectcovs=selectcovs,
		FDR=TRUE
	)

} 
