twobittofa=function(chrm,start,end,twoBitFile,gcSeq){
	.C("twoBitToFa",as.character(chrm),start,end,as.character(twoBitFile),as.character(gcSeq),PACKAGE="zinba")
}
	
