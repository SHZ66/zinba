basealigncount=function(inputfile,outputfile,twoBitFile,extension,filetype="bowtie"){
    cReturn <- .C("baseAlignCounts",as.character(inputfile), as.character(outputfile), as.character(twoBitFile),as.integer(extension),as.character(filetype),PACKAGE="zinba")
}
