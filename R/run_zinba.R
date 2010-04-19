run.zinba=function(filelist=NULL,formula=NULL,formulaE=NULL,outfile=NULL,seq=NULL,align=NULL,input="none",twoBit=NULL,winSize=500,offset=0,cnvWinSize=100000,cnvOffset=0,basecountfile=NULL,threshold=0.01,peakconfidence=.8,tol=10^-5,numProc=1,buildwin=1,pWinSize=200,pquant=0.75,refinepeaks=1,printFullOut=0,method='pscl',initmethod='count',diff=0,filetype="bowtie",extension, cleanup=FALSE){
        rmc <- require(multicore)
        rdmc <- require(doMC)
        rfor <- require(foreach)
        if(rmc == FALSE){
            stop(paste("multicore available"))
        }
	time.start <- Sys.time()
	if(is.null(formulaE)){
		formulaE=formula
	}

        #####################################################################################################
	#create subdirectory to hold intermediate files to be used later
	outfile_subdir=paste(outfile,"_files/", sep="")
	if(!dir.create(outfile_subdir, showWarnings=F)){
		#if cannot create directory, delete existing and try to create again
		cat(paste("\nOverwriting previously existing directory ",outfile_subdir, "\n",sep=""))
		unlink(outfile_subdir, recursive=T)
		if(!dir.create(outfile_subdir,showWarnings=F)){
			#if fails again, print error and set output directory to same as main files at "outpath"
			cat("\nCould not create subdirectory to hold intermediate files, placing intermediate files with main output files\n")
			outfile_subpath=outfile
		}else{
			#if retry successful, set path prefix to where built window analysis files will be place in subdirectory			
			slashindex=which(substring(outfile,1:nchar(outfile),1:nchar(outfile))=="/")
			if(length(slashindex>0)){
				outfile_subpath=paste(outfile,"_files","/",substr(outfile,slashindex[length(slashindex)]+1,nchar(outfile)),sep="")
			}else{
				outfile_subpath=paste(outfile,"_files","/",outfile,sep="")	
			}
		}
	}else{
		#set path prefix to where built window analysis files will be place in subdirectory			
		slashindex=which(substring(outfile,1:nchar(outfile),1:nchar(outfile))=="/")
		if(length(slashindex>0)){
			outfile_subpath=paste(outfile,"_files","/",substr(outfile,slashindex[length(slashindex)]+1,nchar(outfile)),sep="")
		}else{
			outfile_subpath=paste(outfile,"_files","/",outfile,sep="")	
		}
	}
        #####################################################################################################
	if(buildwin==1){
	    if(is.null(filelist)) filelist=paste(outfile_subpath,".list",sep="")	
	    cat(paste("\n--------BEGIN BUILDING WINDOW DATA--------",as.character(Sys.time()),"\n"))
            buildwindowdata(seq=seq,align=align,input=input,twoBit=twoBit,winSize=winSize,offset=offset,cnvWinSize=cnvWinSize,cnvOffset=cnvOffset,filelist=filelist,filetype=filetype,extension=extension, outdir=outfile_subdir)
	}
	if(refinepeaks==1 && is.null(basecountfile)){
		stop(paste("Basecount file must be specified, currently",basecountfile,sep=" "))
	}else if (is.null(filelist)){
		stop(paste("Need list of files ",filelist,sep=" "))
	}else if(method != 'pscl' && method != 'mixture'){
		stop(paste("Method should be either pscl or mixture, currently",method))
	}else{
            params=scan(filelist,what=character(0),quiet=T)
	    winlist=paste(outfile_subpath,".winlist",sep="")
            peakout=paste(outfile,".peaks",sep="")
            bpout=paste(outfile_subpath,".bpcount",sep="")

            if(rmc == TRUE && rdmc == TRUE && rfor == TRUE){
	    cat(paste("--------GETTING ENRICHED WINDOWS--------",as.character(Sys.time()),"\n\n")) 		
                registerDoMC(numProc)
                mcoptions <- list(preschedule = FALSE, set.seed = FALSE)
                getDoParWorkers()
                winfiles <- foreach(i=1:length(params),.combine='rbind',.inorder=FALSE,.errorhandling="remove",.options.multicore = mcoptions) %dopar%
                    getsigwindows(file=params[i],formula=formula,formulaE=formulaE,threshold=threshold,winout=outfile_subpath,peakconfidence=peakconfidence,tol=tol,method=method,printFullOut=printFullOut,initmethod=initmethod)
            
                write.table(winfiles,winlist,quote=F,row.names=F,col.names=F)
	    cat(paste("--------WINDOW ANALYSIS COMPLETE--------",as.character(Sys.time()),"\n\n"))		
            }else{
	    cat(paste("--------GETTING ENRICHED WINDOWS--------",as.character(Sys.time()),"\n\n")) 	
                for(i in 1:length(params)){
                    wfile <- getsigwindows(file=params[i],formula=formula,formulaE=formulaE,threshold=threshold,winout=outfile,peakconfidence=peakconfidence,tol=tol,method=method,printFullOut=printFullOut,initmethod=initmethod)
                    winfiles <- rbind(wfile)
                }
	    cat(paste("--------WINDOW ANALYSIS COMPLETE--------",as.character(Sys.time()),"\n\n"))		
            }
	    if(refinepeaks==1){
		cat(paste("--------MERGE WINDOWS AND REFINE PEAKS--------",as.character(Sys.time()),"\n"))
		getrefinedpeaks(winlist=winlist,basecountfile=basecountfile,bpout=bpout,peakout=peakout,twoBit=twoBit,winSize=winSize,pWinSize=pWinSize,pquant=pquant,printFullOut=printFullOut,peakconfidence=peakconfidence,threshold=threshold,method=method)
	    }
		
	}
        #####################################################################################################
	cat(paste("\n--------ZINBA COMPLETE--------",as.character(Sys.time()),"\n\n"))
	if(cleanup==TRUE) unlink(outfile_subdir, recursive=T)
	time.end <- Sys.time()
	print(difftime(time.end,time.start))
}
