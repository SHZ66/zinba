buildwindowdata=function(seq, align, gdna, cnvarray){
	#need to update perl file by removing environment for hg18 path, incorporate into perl command
	Perl.Path <- file.path(.path.package("zimba"), "Perl")
	Fn.Path <- file.path(Perl.Path, "buildPeakWindows_v3.pl")
        CMD=paste("--seq", seq,"--align", align,"--gdna", gdna,"--cnvarray", cnvarray,sep=" ")
	system(CMD)    
}
