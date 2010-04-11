signalnoise=function(inputfile,output,twoBitFile,binSize=2000){
	c=.C("signalnoise",as.character(inputfile),as.character(output),as.character(twoBitFile),as.integer(binSize),PACKAGE="zinba")
	a=read.table(output)
	median=a[,1]
	max=a[,2]	
	print("Summary Information of log ratio of Max Window Count vs Median Window Count")
	summary(log(max/median)[median>0])
}
	
