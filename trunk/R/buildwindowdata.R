buildwindowdata=function(seq, align=NULL, input=NULL, cnvarray=NULL,twoBit=NULL,winSize=500,offset=0,aThresh=1,gb=NULL){
	Perl.Path <- file.path(.path.package("zinba"), "exec")
	Fn.Path <- file.path(Perl.Path, "buildPeakWindows_v3.pl")
        CMD=paste(Fn.Path,"--seq", seq,"--align", align,"--cnvarray", cnvarray,"--twoBit",twoBit,"--offset-size",offset,"--align-thresh",aThresh,"--gb",gb,sep=" ")
#        CMD=paste(Fn.Path,"--seq", seq,"--align", align,"--input", input,"--cnvarray", cnvarray,"--twoBit",twoBit,"--window-size",winSize,"--offset-size",offset,"--align-thresh",aThresh,"--gb",gb,sep=" ")
	system(CMD)    
}