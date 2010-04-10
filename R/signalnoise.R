signalnoise=function(basecountfile, outfile, winsize){
	.C("signalnoise",as.character(basecountfile), as.character(outfile), as.integer(winsize),PACKAGE="zinba")
}
	
