alignAdjust=function(inputfile, outDir="",twoBitFile,athresh, adjustsize){
    .C("alignAdjust",as.character(inputfile),as.character(outDir),as.character(twoBitFile),as.integer(athresh),as.integer(adjustsize),PACKAGE="zinba")
    gc()
}
