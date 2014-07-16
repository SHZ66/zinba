run = function(y, X, XE = NULL, XZ = NULL, prop1,glmtype="NB", EMmaxit=50, model){
  if(model=="FMR"){
    res = mr(y=y, X=XB, XE = XE, XZ = XZ, prop1 = bgprops[i],glmtype="NB", EMmaxit = 5)
  }else if(model == "HMM"){
    res = hmm(y=y, X=XB, XE = XE, prop1 = bgprops[i],glmtype="NB", EMmaxit = 5)
  }else if(model =="ARHMM"){
    res = arhmm(y=y, X=XB, XE = XE, prop1 = bgprops[i],glmtype="NB", EMmaxit = 5)
  }else{
    stop(sprintf("Model must be FMR, HMM, or ARHMM, specified as %s", model))
  }
  return(res)
}

getsigwindows_quickrun=function(file,formulaE,formulaZ,winout,
                       threshold=.01,peakconfidence=.8, printFullOut=0,tol=10^-5,
                       method="mixture",initmethod="count", diff=0,modelselect=FALSE, trace=0, 
                       FDR=FALSE, model = "MR")
{
  
  time.start <- Sys.time()
  library(MASS)
  library(hmmcov)
  
  options(scipen=999)
  
  
  #quick run
  #adjust enrichment input
  if(eadjust==T){
    formulaE = formulaB = exp_count~input_count
    formulaZ = exp_count~1
  }else{
    formulaE = NULL
    formulaB = exp_count~input_count
    formulaZ = exp_count~1    
  }
  
  
  #get files
  files = list.files(projdir,".txt", full.names=T)
  
  while(fnum <= length(files)){
      # read in data
      data=read.table(files[fnum], header=TRUE)
      
      # set chr
      chrm = data$chromosome[1]
      
      # build covariate matrices
      mf <- model.frame(formula=formula, data=data)
      mfZ <- model.frame(formula=formulaZ, data=data)
      X <- model.matrix(attr(mf, "terms"), data=mf)
      XZ <- model.matrix(attr(mfZ, "terms"), data=mfZ)
      
      if(eadjust==T){
        mfE <- model.frame(formula=formulaE, data=data)
        XE <- model.matrix(attr(mfE, "terms"), data=mfE)
      }else{
        XE = NULL
      }
      
      #get response data
      Y <- model.response(mf)
      
      #inialization of model, choosing prop1
      if(fnum == 1){
        bgprops = c(.99, .95, .9, .8, .7)
        ll = rep(length(bgprops))
        parms = list()
        for(i in 1:length(bgprops)){
          res = run(y=y, X=XB, XE = XE, XZ = XZ, prop1 = bgprops[i],glmtype="NB", EMmaxit = 5, model = model)
          ll[i] = res$ll
          parms[[i]] = res$coef
        }
        prop = bgprops[which.max(ll)]
        #recreta pp
      }
      
      #run model, use above to initialize
      res = run(y=y, X=XB, XE = XE, XZ = XZ, prop1 = prop, glmtype="NB", model = model, Pi, forwardback)
      
      #get posterior probs, convert to FDR
      p  = 1-res$forwardback[,2]
      p2 = rep(0,length(p))
      p2[order(p)] = cumsum(p[order(p)])/(1:length(p))
      
      numpeaks = length(which(p2<threshold))
      
      line = paste("\nProcessed ",files[fnum],"\n\t","Selected number of peaks: ",as.character(numpeaks),"\n\t",as.character(Sys.time()),"\t",sep='')
      
      winfile = paste(winout,"_",chrm,".wins",sep="")
      write.table(data,winfile,quote=F,sep="\t",row.names=F)
      fnum=fnum+1
      }
    
  
  time.end <- Sys.time()
  cat("\n")
  print(difftime(time.end,time.start))
  return(winfile)
  #final
}
