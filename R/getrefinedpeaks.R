getrefinedpeaks=function(winlist,coordout,basecountfile,bpout,peakout,twoBit,winSize,winGap=4,pWinSize=200,printFullOut=0,pquant=0.75,peakconfidence=.8,threshold=.01,method='pscl'){
    if(method=='pscl'){
        createcoordfile(filelist=winlist,threshold=threshold,coordout=coordout,method=method,winSize=winSize,winGap=winGap,printFullOut=printFullOut)
    }else if(method=='mixture'){
        createcoordfile(filelist=winlist,threshold=peakconfidence,coordout=coordout,method=method,winSize=winSize,winGap=winGap,printFullOut=printFullOut)
    }
    basecountimport(inputfile=basecountfile,coordfile=coordout,outputfile=bpout,twobitfile=twoBit)
    peakbound(bpprofile=bpout,output=peakout,pwinSize=pWinSize,winSize=winSize,quantile=pquant)
    unlink(coordout)
    unlink(bpout)
}
