estimatecnv=function(seq,gb=NULL,winSize=25000,offset=2500,twoBit=NULL,nProcess=1){
        library(MASS)
        options(scipen=999)

	Perl.Path <- file.path(.path.package("zinba"), "exec")
	Fn.Path <- file.path(Perl.Path, "buildPeakWindows_v3.pl")
        CMD=paste(Fn.Path,"--seq", seq,"--align", align,"--twoBit",twoBit,"--window-size",winSize,"--offset-size",offset,"--perc-n-thresh",percN,"--gb",gb,"--processes",nProcess,sep=" ")
	system(CMD)
	cnv_file <- paste(substr(seq,1,((nchar(as.character(seq)))-4)),".cnvs",sep="")

	cnv_data=read.table(cnv_file,header=F)
	winout <- paste(substr(file_list[i,1],1,((nchar(as.character(file_list[i,1])))-4)),".cnvs",sep="")
	data=read.table(as.character(file_list[i,1]),header=TRUE,sep="\t")
	colnames(data)=c('chromsome','start','stop','exp_count')
        
	a=glm.nb(data$exp_count ~ 1)        
	cnv_vals=cbind(data,residuals(a))
	colnames(cnv_vals)[dim(cnv_vals)[2]]=c('residuals')
	write.table(cnv_vals,winout,quote=F,sep="\t",row.names=F,col.names=T)
}