getsigwindows2=function(file,formula,formulaE,formulaZ,winout,
                       threshold=.01,peakconfidence=.8, printFullOut=0,tol=10^-5,
                       method="mixture",initmethod="count", diff=0,modelselect=FALSE, trace=0, 
                       FDR=FALSE)
{
  
  time.start <- Sys.time()
  library(MASS)
  rhmmcov = require(hmmcov)
  if(rhmmcov == F) stop("hmmcov package is needed")
  options(scipen=999)
  
  if(!inherits(formula, "formula"))  
    stop("Check your background component formula, not entered as a formula object")
  if(!inherits(formulaE, "formula")) 
    stop("Check your enrichment component formula, not entered as a formula object")
  if(!inherits(formulaZ, "formula")) 
    stop("Check your zero-inflated component formula, not entered as a formula object")
  
  winfile=NULL
  printflag = 0
  files = unlist(strsplit(file,";"))
  fnum=1
  retry=0	#flag for retrial of offset
    while(fnum <= length(files)){
      fail=0 #failure flag for GLM, reset to zero for each offset
      #print(files[fnum])
      data=read.table(files[fnum], header=TRUE)
      chrm = data$chromosome[1]
      q25=quantile(data$exp_count, 0) #not used, set to 0
      mf <- model.frame(formula=formula, data=data)
      mfE <- model.frame(formula=formulaE, data=data)
      mfZ <- model.frame(formula=formulaZ, data=data)
      X <- model.matrix(attr(mf, "terms"), data=mf)
      XE <- model.matrix(attr(mfE, "terms"), data=mfE)
      XZ <- model.matrix(attr(mfZ, "terms"), data=mfZ)
      
      Y <- model.response(mf)
      
      #initalization
      if(fnum == 1){ 
        props = c(.999, .99, .95, .9, .8, .7)
        ll0  = rep(-Inf, length(props))
        for(k in 1:length(props)){
          ll0[k] = fmr(y=Y, X=as.matrix(X[,-1]), prop1=props[k], XE=as.matrix(XE[,-1]),maxitEM=5, glmtype="nb", zeroinfl = T, XZ = as.matrix(XZ[,-1]), trace = trace, thresh = 0)$ll
          if(trace==1) cat("Initialization prop ", k, "\n")
        }
        prop = props[which.max(ll0)]
      }
      if(trace==1) cat("File ", files[fnum], "\n")
      
      result = fmr(y=Y, X=as.matrix(X[,-1]), prop1=prop, XE=as.matrix(XE[,-1]),maxitEM=100, glmtype="nb", zeroinfl = T, XZ = as.matrix(XZ[,-1]), trace = trace)

        probi2 = result$forwardbackward[,2]
        p=1-probi2
        p2=rep(0,length(p))
        p2[order(p)]=cumsum(p[order(p)])/(1:length(p))
        if(FDR == T){
          numpeaks=length(which(p2<threshold))
        }else{
          numpeaks=length(which(probi2>peakconfidence))
        }
        
        if(printFullOut == 1){
          data=cbind(data,((data$exp_count>q25)^2),formatC(probi2,format="f",digits=16), formatC(p2,format="f",digits=16)) #all non-zero set to 1
          colnames(data)[c(dim(data)[2]-2 ,dim(data)[2]-1,dim(data)[2])]=c('q25','peakprob','qvalue')
        }else{
          data=cbind(data[1:3],((data$exp_count>q25)^2),formatC(probi2,format="f",digits=16), formatC(p2,format="f",digits=16))
          colnames(data)=c('chromosome','start','stop','q25','peakprob','qvalue') #all non-zero set to 1
        }
        #if(diff==0){
          #if differential expression (comparing two samples) is not occuring
        #  data$q25[Y-mui1<0]=0
        #}
        line = paste("\nProcessed ",files[fnum],"\n\t","Selected number of peaks: ",as.character(numpeaks),"\n\t",as.character(Sys.time()),"\t",sep='')
        ### PRINT WINDOWS
        cat(line)
        winfile = paste(winout,"_",chrm,".wins",sep="")
        if(printflag==1){
          write.table(data,winfile,quote=F,sep="\t",row.names=F,col.names=F,append=T)
        }else{
          write.table(data,winfile,quote=F,sep="\t",row.names=F)
          printflag=1
        }
        fnum=fnum+1
        rm(data); rm(Y); rm(X); rm(XE); rm(XZ); rm(probi2); gc();
    }
  
  time.end <- Sys.time()
  cat("\n")
  print(difftime(time.end,time.start))
  return(winfile)
  #final
}
