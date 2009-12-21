getrefinedpeaks=function(winlist,basecountfile,bpout,peakout,twoBit,winSize,pWinSize=200,printFullOut=0,pquant=0.75,threshold=.01,peakconfidence=0.8,method='pscl'){
    if(method=='mixture'){
        threshold = peakconfidence
    }
    basecountimport(inputfile=basecountfile,winlist=winlist,threshold=threshold,method=method,printFullOut=printFullOut,outputfile=bpout,twobitfile=twoBit)
    peakbound(bpprofile=bpout,output=peakout,pwinSize=pWinSize,winSize=winSize,quantile=pquant)
    unlink(bpout)
}
