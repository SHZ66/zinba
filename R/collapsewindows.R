collapsewindows=function(winlist,thresholds,method='pscl',printFullOut=0){
    cReturn <- .C("collapse_windows",as.character(winlist),as.character(method),as.integer(printFullOut),as.integer(length(thresholds)),as.double(thresholds),PACKAGE="zinba")
}
