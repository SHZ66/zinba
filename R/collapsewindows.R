collapsewindows=function(winlist,thresholds,method='pscl',printFullOut=0){
	if(!file.exists(mapdir)){stop("Specified winlist doesnt exist or is incorrect")}
	cReturn <- .C("collapse_windows",as.character(winlist),as.character(method),as.integer(printFullOut),as.integer(length(thresholds)),as.double(thresholds),PACKAGE="zinba")
}
