buildwindowdata=function(seq=NULL, align=NULL, input=NULL, cnvarray=NULL,cnvexp=NULL,rand=NULL,twoBit=NULL,winSize=500,offset=0,aThresh=1,inputTrans=NULL,percN=0.1,gb="hg18",nProcess=1){
	Perl.Path <- file.path(.path.package("zinba"), "exec")
	Fn.Path <- file.path(Perl.Path, "buildPeakWindows_v3.pl")
        CMD=paste(Fn.Path,"--seq", seq,sep=" ")
	if(is.null(seq)){
		print(paste("Seq must be defined, currently:",seq,sep=" "))
	}else{
		if (!is.null(align)){
			CMD=paste(CMD,"--align",align,sep=" ")
		}
		if (!is.null(input)){
			CMD=paste(CMD,"--input",input,sep=" ")
		}
		if (!is.null(cnvarray)){
			CMD=paste(CMD,"--cnvarray",cnvarray,sep=" ")
		}
		if (!is.null(cnvexp)){
			CMD=paste(CMD,"--cnv-exp",cnvexp,sep=" ")
		}
		if (!is.null(rand)){
			CMD=paste(CMD,"--rand",rand,sep=" ")
		}
		if (!is.null(twoBit)){
			CMD=paste(CMD,"--twoBit",twoBit,sep=" ")
		}
		CMD=paste(CMD,"--window-size",winSize,"--offset-size",offset,"--align-thresh",aThresh,"--perc-n-thresh",percN,"--gb",gb,"--processes",nProcess,sep=" ")
		if (!is.null(inputTrans)){
			CMD=paste(CMD,"--trans-input",sep=" ")
		}
		system(CMD)
	}
}