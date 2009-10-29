basealigncount=function(inputfile, outputfile, extension){
    cReturn <- .C("baseAlignCounts",as.character(inputfile), as.character(outputfile), as.integer(extension),PACKAGE="zinba")
}
