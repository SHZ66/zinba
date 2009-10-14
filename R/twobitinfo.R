twobitinfo=function(infile, outfile){
        print(paste("infile is",infile,"outfile is",outfile,sep=" "))
	.C("twoBitInfo",as.character(infile), as.character(outfile),PACKAGE="zinba")
}
	
