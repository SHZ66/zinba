twobittofa=function(chrm,start,end,twoBitFile,gcSeq){
    .C("twoBitToFa",as.character(chrm),as.integer(start),as.integer(end),as.character(twoBitFile),as.character(gcSeq),PACKAGE="zinba")
}
	
