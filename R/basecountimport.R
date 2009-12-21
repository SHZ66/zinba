basecountimport=function(inputfile,winlist,threshold=.01,method='pscl',printFullOut=0,outputfile,twobitfile,chromosome='all'){
    cReturn <- .C("getSeqCountProfile",as.character(inputfile),as.character(winlist),as.double(threshold),as.character(method),as.integer(printFullOut),as.character(outputfile),as.character(twobitfile),as.character(chromosome),PACKAGE="zinba")
}
