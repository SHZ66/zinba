getrefinedpeaks=function(winlist,coordout,basecountfile,bpout,peakout,twoBit,pWinSize,pquant=0.75,peakconfidence=.8,threshold=.01,method='pscl'){
    createcoordfile(filelist=winlist,threshold=threshold,peakconfidence=peakconfidence,coordout=coordout,method=method)
    basecountimport(inputfile=basecountfile,coordfile=coordout,outputfile=bpout,twobitfile=twoBit)
    peakbound(bpprofile=bpout,output=peakout,winSize=pWinSize,quantile=pquant)
    unlink(coordout)
    unlink(bpout)
}