twobittofa=function(inName, outName){
	.C("twoBitToFa.c",as.character(inName), as.character(outName),PACKAGE="zinba")
}
	
