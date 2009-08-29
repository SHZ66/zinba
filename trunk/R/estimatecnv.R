estimatecnv=function(seq,winSize=25000,offset=0,twoBit=NULL,nProcess=1,percN=0.1){
        library(MASS)
        options(scipen=999)

	Perl.Path <- file.path(.path.package("zinba"), "exec")
	Fn.Path <- file.path(Perl.Path, "buildCNVWindows.pl")
        CMD=paste(Fn.Path,"--seq", seq,"--twoBit",twoBit,"--window-size",winSize,"--offset-size",offset,"--perc-n-thresh",percN,"--processes",nProcess,sep=" ")
	system(CMD)
	cnv_file <- paste(substr(seq,1,((nchar(as.character(seq)))-5)),".cnvs",sep="")

	cnv_data=read.table(cnv_file,header=F)
	winout <- paste(substr(seq,1,((nchar(as.character(seq)))-5)),"_RESIDUAL.cnvs",sep="")
	data=read.table(as.character(cnv_file),header=TRUE,sep="\t")
	colnames(data)=c('chromsome','start','stop','exp_count')
        
	a=glm.nb(data$exp_count ~ 1)        
	cnv_vals=cbind(data,residuals(a))
	colnames(cnv_vals)[dim(cnv_vals)[2]]=c('residuals')
	write.table(cnv_vals,winout,quote=F,sep="\t",row.names=F,col.names=T)
}