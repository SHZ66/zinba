getrefinedpeaks=function(winout,coordout,basecountfile,bpout,peakout,twoBit,pWinSize,pquant,peakconfidence=.8,threshold=.01,method='pscl'){
    createcoordfile(file=winout,threshold=threshold,peakconfidence=peakconfidence,coordout=coordout,method=method)
print(paste("Getting base count data"))
    basecountimport(inputfile=basecountfile,coordfile=coordout,outputfile=bpout,twobitfile=twoBit)
print(paste("Getting peakbounds"))
    peakbound(bpprofile=bpout,output=peakout,winSize=pWinSize,quantile=pquant)
print(paste("Finished"))
    #unlink(coordout)
    #unlink(bpout)
}