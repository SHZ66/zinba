crosscorrelation=function(inputfile,twoBitFile,filetype){
    cReturn <- .C("crossCor",as.character(inputfile),as.character(twoBitFile),as.character(filetype),PACKAGE="zinba")
}
