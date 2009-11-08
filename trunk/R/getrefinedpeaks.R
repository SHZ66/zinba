getrefinedpeaks=function(winout,coordout,basecountfile,bpout,peakout,twoBit,pWinSize,pquant,peakconfidence=.8,threshold=.01,method='pscl'){
    createcoordfile(file=winout,threshold=threshold,peakconfidence=peakconfidence,coordout=coordout,method=method)
    basecountimport(inputfile=basecountfile,coordfile=coordout,outputfile=bpout,twobitfile=twoBit)
    peakbound(bpprofile=bpout,output=peakout,winSize=pWinSize,quantile=pquant)
    unlink(coordout)
    unlink(bpout)
}