alignAdjust=function(inputfilelist, athresh, adjustsize){
    .C("alignAdjust",as.character(inputfilelist), as.integer(athresh), as.integer(adjustsize),PACKAGE="zinba")
}
