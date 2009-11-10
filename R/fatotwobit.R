fatotwobit=function(faFile,outFile){
    .C("faToTwoBit",as.character(faFile),as.character(outFile),PACKAGE="zinba")
}
	
