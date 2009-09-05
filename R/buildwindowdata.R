buildwindowdata=function(seq, align="NULL", input="NULL", cnvarray="NULL",cnvexp="NULL",rand="NULL",twoBit="NULL",winSize=500,offset=0,aThresh=1,inputTrans=NA,percN=0.1,gb="NULL",nProcess=1){
	Perl.Path <- file.path(.path.package("zinba"), "exec")
	Fn.Path <- file.path(Perl.Path, "buildPeakWindows_v3.pl")
	transIn=NULL
	if (!is.na(inputTrans)){
		transIn = "--trans-input"
	}
        CMD=paste(Fn.Path,"--seq", seq,"--align", align,"--cnvarray", cnvarray,"--cnv-exp",cnvexp,"--rand",rand,"--twoBit",twoBit,"--window-size",winSize,"--offset-size",offset,"--align-thresh",aThresh,"--perc-n-thresh",percN,"--gb",gb,"--processes",nProcess,transIn,sep=" ")
	system(CMD)
}