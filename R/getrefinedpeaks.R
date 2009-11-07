getrefinedpeaks=function(winout,coordout,basecountfile,bpout,twoBit,pWinSize,pquant,threshold=.01,method='pscl'){
    createcoordfile(file=winout,threshold=threshold,coordout=coordout,method=method)
    basecountimport(inputfile=basecountfile,coordfile=coordout,outputfile=bpout,twobitfile=twoBit)
    peakbound(bpprofile=bpout,output=peakout,winSize=pWinSize,quantile=pquant)
    #unlink(coordout)
    #unlink(bpout)
}