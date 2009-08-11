buildgcdata=function(listFA,gb){
	Perl.Path <- file.path(.path.package("zinba"), "src")
	Fn.Path <- file.path(Perl.Path, "faTobinary.pl")
        CMD=paste("Fn.Path",listFA,gb,sep=" ")
	system(CMD)
}
