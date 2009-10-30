basealigncount=function(inputfile, outputfile, twoBitFile,extension){
    cReturn <- .C("baseAlignCounts",as.character(inputfile), as.character(outputfile), as.character(twoBitFile),as.integer(extension),PACKAGE="zinba")
}
