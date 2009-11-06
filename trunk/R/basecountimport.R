basecountimport=function(inputfile, coordfile, outputfile, twobitfile, chromosome='all'){
    cReturn <- .C("getSeqCountProfile",as.character(inputfile), as.character(coordfile), as.character(outputfile), as.character(twobitfile),as.character(chromosome),PACKAGE="zinba")
}
