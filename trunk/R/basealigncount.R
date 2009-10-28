basealigncount=function(inputfile, outputfile, extension){
    cReturn <- .C("baseAlignCounts",as.character(inputfile), as.character(coordfile), as.integer(extension),PACKAGE="zinba")
}
