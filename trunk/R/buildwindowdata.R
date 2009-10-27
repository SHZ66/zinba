buildwindowdata=function(seq,input="none",align,twoBit,winSize=500,offset=0,cnvWinSize=100000,cnvOffset=0,paramFile){
	cReturn <- .C("buildWindows",as.character(seq),as.character(input),as.character(align),as.character(twoBit),as.integer(winSize),as.integer(offset),as.integer(cnvWinSize),as.integer(cnvOffset),as.character(paramFile),PACKAGE="zinba")
}
#buildwindowdata=function(paramFile=NULL,seq=NULL, align=NULL, input=NULL, cnvarray=NULL,cnvexpwin=NULL,cnvexpcustom=NULL,rand=NULL,twoBit=NULL,winSize=500,offset=0,inputTrans=NULL,inputRandLog=NULL,percN=0.1,gb="hg18",nProcess=1){
#	Perl.Path <- file.path(.path.package("zinba"), "exec")
#	Fn.Path <- file.path(Perl.Path, "buildPeakWindows_v3.pl")
#        CMD=paste(Fn.Path,"--seq", seq,sep=" ")
#	if(is.null(seq)){
#		print(paste("Seq must be defined, currently:",seq,sep=" "))
#	}else{
#		if (!is.null(align)){
#			CMD=paste(CMD,"--align",align,sep=" ")
#		}
#		if (!is.null(input)){
#			CMD=paste(CMD,"--input",input,sep=" ")
#		}
#		if (!is.null(cnvarray)){
#			CMD=paste(CMD,"--cnvarray",cnvarray,sep=" ")
#		}
#		if (!is.null(cnvexpwin)){
#			CMD=paste(CMD,"--cnv-expwin",cnvexpwin,sep=" ")
#		}
#		if (!is.null(cnvexpcustom)){
#			CMD=paste(CMD,"--cnv-expcustom",cnvexpcustom,sep=" ")
#		}
#		if (!is.null(rand)){
#			CMD=paste(CMD,"--rand",rand,sep=" ")
#		}
#		if (!is.null(twoBit)){
#			CMD=paste(CMD,"--twoBit",twoBit,sep=" ")
#		}
#		if (!is.null(paramFile)){
#			CMD=paste(CMD,"--param-file",paramFile,sep=" ")
#		}
#		CMD=paste(CMD,"--window-size",winSize,"--offset-size",offset,"--perc-n-thresh",percN,"--gb",gb,"--processes",nProcess,sep=" ")
#		if (!is.null(inputTrans)){
#			CMD=paste(CMD,"--trans-input",sep=" ")
#		}
#		if (!is.null(inputRandLog)){
#			CMD=paste(CMD,"--input_rand_log2",sep=" ")
#		}
#		system(CMD)
#	}
#}