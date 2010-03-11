getrefinedpeaks=function(winlist,basecountfile,bpout,peakout,twoBit,winSize,pWinSize=200,printFullOut=0,pquant=0.75,threshold=.01,peakconfidence=0.8,method='pscl',minscore=0){
    if(method=='mixture'){
        threshold = peakconfidence
    }
    print(paste("Threshold is ",threshold,sep=""))
    basecountimport(inputfile=basecountfile,winlist=winlist,threshold=threshold,method=method,printFullOut=printFullOut,outputfile=bpout,twobitfile=twoBit)
    peakbound(bpprofile=bpout,output=peakout,pwinSize=pWinSize,winSize=winSize,quantile=pquant,minscore=minscore)
    unlink(bpout)
}