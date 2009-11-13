createcoordfile=function(filelist,threshold=.01,peakconfidence=.8,coordout,method='pscl'){
    files=scan(filelist,what=character(0))
    for(f in 1:length(files)){
        data=read.table(files[f],header=T)
        numpeaks=NULL
        sigpeaks=NULL
        if(method=='pscl'){
            numpeaks=length(which(data$qvalue<threshold))
            sigpeaks=data[which(data$qvalue<threshold),]
            print(paste(as.character(numpeaks),' windows were less than ',as.character(threshold), sep=''))
        }else{
            numpeaks=length(which(data$peakprob>peakconfidence))
            sigpeaks=data[which(data$peakprob>peakconfidence),]
            print(paste(as.character(numpeaks),' windows were greater than ',as.character(peakconfidence), sep=''))
        }

        peakID=paste(sigpeaks$chromosome,sigpeaks$start,sigpeaks$stop,sep=":")
        coordinates=cbind(peakID,as.character(sigpeaks$chromosome),sigpeaks$start,sigpeaks$stop,sigpeaks$q25,"+")
        if(file.exists(coordout)){
            write.table(coordinates,coordout,quote=F,append=TRUE,sep="\t",row.names=F,col.names=F)
        }else{
            write.table(coordinates,coordout,quote=F,sep="\t",row.names=F,col.names=F)
        }
    }
}