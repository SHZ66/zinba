printsimdata=function(output,counts, winsize=500, ext=200, chrom='chr22'){
	reads=sum(counts[,2])
	readfile=matrix(0, reads, 2)
		#non enriched sites
	i=1
	if(counts[i,1]==0){
		readcenter=sort(round(runif(count2[i,2],(i-1)*winsize+1, i*winsize)))
		readfile[1:(sum(counts[i,2])), ]=c(readcenter-(ext/2-1), readcenter+(ext/2-1))
	}else{
		readcenter=round(runif(count2[i,2],(i-1)*winsize+1, i*winsize))
		readfile[1:(sum(counts[i,2])), ]=c(readcenter-(ext/2), readcenter+(ext/2))
	}


	for(i in 2:length(counts[,1])){
		if(counts[i,2]==0){
			i=i+1
		}else{
			readcenter=sort(round(runif(count2[i,2],(i-1)*winsize+1, i*winsize)))
			readfile[(sum(counts[1:(i-1),2])+1):(sum(counts[1:(i-1),2])+counts[i,2]), ]=c(readcenter-(ext/2-1), 				readcenter+(ext/2-1))
			i=i+1
		}
	}
	readfile[which(readfile[,1]<0),1]=1
	readfile[which(readfile[,2]>winsize*dim(counts)[1]),2]=winsize*dim(counts)[1]
	write.table(list(rep(chrm, dim(readfile)[1]),readfile),output, sep='\t',row.names = F,col.names = F)
}



