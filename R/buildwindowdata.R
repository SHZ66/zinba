buildwindowdata=function(seq, align=NULL, input=NULL, cnvarray=NULL,twoBit=NULL,winSize=500,offset=0,aThresh=1,percN=0.1,gb=NULL){
	Perl.Path <- file.path(.path.package("zinba"), "exec")
	Fn.Path <- file.path(Perl.Path, "buildPeakWindows_v3.pl")
        CMD=paste(Fn.Path,"--seq", seq,"--align", align,"--cnvarray", cnvarray,"--twoBit",twoBit,"--win-size",winSize,"--offset-size",offset,"--align-thresh",aThresh,"--perc-n-thresh",percN,"--gb",gb,sep=" ")
	system(CMD)    
}