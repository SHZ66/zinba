estimatecnv=function(seq=NULL,winFile=NULL,winSize=25000,offset=0,gb="hg18",nProcess=1,alignDir=NULL){
        options(scipen=999)

	Perl.Path <- file.path(.path.package("zinba"), "exec")
	Fn.Path <- file.path(Perl.Path, "buildCNVWindows.pl")

        CMD=paste(Fn.Path,"--seq", seq,"--window-size",winSize,"--offset-size",offset,"--gb",gb,"--processes",nProcess,sep=" ")
	if(is.null(seq)){
		print(paste("Seq must be defined, currently:",seq,sep=" "))
	}else{
		if (!is.null(winFile)){
			CMD=paste(CMD,"--window_file",winFile,sep=" ")
		}
		if (!is.null(alignDir)){
			CMD=paste(CMD,"--align-dir",alignDir,sep=" ")
		}
		system(CMD)
	}
}