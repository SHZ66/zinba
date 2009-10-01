basecountimport=function(inputfile, coordfile, outputfile, chromosome='all'){
    .C("getSeqCountProfile",as.character(inputfile), as.character(coordfile), as.character(outputfile), as.character(chromosome),PACKAGE="zinba")
    gc()
}
