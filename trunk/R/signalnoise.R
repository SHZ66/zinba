signalnoise=function(inputfile,output,twoBitFile,binSize=500){
	c=.C("signalnoise",as.character(inputfile),as.character(output),as.character(twoBitFile),as.integer(binSize),PACKAGE="zinba")
}
	
