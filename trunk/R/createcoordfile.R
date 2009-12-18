#createcoordfile=function(filelist,threshold=.01,peakconfidence=.8,coordout,method='pscl'){
createcoordfile=function(filelist,threshold=.01,coordout,method='pscl',winSize,winGap=4,printFullOut=0){
    cReturn <- .C("mergeWins",as.character(filelist),as.character(coordout),as.integer(winSize),as.integer(winGap),as.double(threshold),as.integer(printFullOut),as.character(method),PACKAGE="zinba")

    #files=scan(filelist,what=character(0))
    #numSigWins = 0
    #for(f in 1:length(files)){
    #    if(file.exists(files[f])){
    #        print(paste("Processing: ",files[f]))
    #        data=read.table(files[f],header=T)
    #        numpeaks=NULL
    #        sigpeaks=NULL
    #        if(method=='pscl'){
    #            numpeaks=length(which(data$qvalue<threshold))
    #            sigpeaks=data[which(data$qvalue<threshold),]
    #            numSigWins = numpeaks + numSigWins
    #        }else{
    #            numpeaks=length(which(data$peakprob>peakconfidence))
    #            sigpeaks=data[which(data$peakprob>peakconfidence),]
    #            numSigWins = numpeaks + numSigWins
    #        }
    #
    #        peakID=paste(sigpeaks$chromosome,sigpeaks$start,sigpeaks$stop,sep=":")
    #        coordinates=cbind(peakID,as.character(sigpeaks$chromosome),sigpeaks$start,sigpeaks$stop,sigpeaks$q25,"+")
    #        if(file.exists(coordout)){
    #            write.table(coordinates,coordout,quote=F,append=TRUE,sep="\t",row.names=F,col.names=F)
    #        }else{
    #            write.table(coordinates,coordout,quote=F,sep="\t",row.names=F,col.names=F)
    #        }
    #    }else{
    #        print(paste("Following file does not exist:",files[f]))
    #    }
    #}
    #print(paste(as.character(numSigWins),' windows were greater than ',as.character(peakconfidence), sep=''))
}