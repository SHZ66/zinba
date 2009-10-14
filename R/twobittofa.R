twobittofa=function(chrm,start,end,twoBitFile,gcSeq){
    print(paste("chrm start end twoBitFile gcSeq are",chrm,start,end,twoBitFile,gcSeq,sep=" "))
    .C("twoBitToFa",as.character(chrm),as.integer(start),as.integer(end),as.character(twoBitFile),as.character(gcSeq),PACKAGE="zinba")
}
	
