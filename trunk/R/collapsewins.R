collapsewins=function(winlist,winout){
    winfiles=scan(winlist,what=character(0))
    for(f in 1:length(winfiles)){
        if(file.exists(winfiles[f])){
            data=read.table(winfiles[f],header=T)
            if(file.exists(winout)){
                write.table(data,winout,quote=F,append=TRUE,sep="\t",row.names=F,col.names=F)
            }else{
                write.table(data,winout,quote=F,sep="\t",row.names=F)
            }
            unlink(winfiles[f])
        }
    }
    gc()
}