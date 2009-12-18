createcoordfile=function(filelist,threshold=.01,coordout,method='pscl',winSize,winGap=4,printFullOut=0){
    cReturn <- .C("mergeWins",as.character(filelist),as.character(coordout),as.integer(winSize),as.integer(winGap),as.double(threshold),as.integer(printFullOut),as.character(method),PACKAGE="zinba")
}