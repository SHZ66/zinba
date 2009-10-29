buildwindowdata=function(seq,input="none",align,twoBit,winSize=500,offset=0,cnvWinSize=100000,cnvOffset=0,paramFile){
	cReturn <- .C("buildWindows",as.character(seq),as.character(input),as.character(align),as.character(twoBit),as.integer(winSize),as.integer(offset),as.integer(cnvWinSize),as.integer(cnvOffset),as.character(paramFile),PACKAGE="zinba")
}
