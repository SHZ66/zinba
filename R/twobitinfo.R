twobitinfo=function(infile, outfile){
	.C("twoBitInfo.c",as.character(infile), as.character(outfile),PACKAGE="zinba")
}
	
