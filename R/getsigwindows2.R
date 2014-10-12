run = function(y, X, XE = NULL, XZ = NULL, prop1,glmtype="NB", maxitEM=50, model, diff = F){
  if(model=="FMR" & diff == F){
    res = fmr(y=y, X=XB, XE = XE, XZ = XZ, prop1 =prop1,glmtype="NB", maxitEM = maxitEM, diff = diff )
  }else if(model == "HMM"){
    res = hmm(y=y, X=XB, XE = XE, prop1 = prop1,glmtype="NB", maxitEM = maxitEM, diff= diff )
  }else if(model =="ARHMM"){
    res = arhmm(y=y, X=XB, XE = XE, prop1 = prop1,glmtype="NB", maxitEM = maxitEM, diff= diff )
  }else{
    stop(sprintf("Model must be FMR, HMM, or ARHMM, specified as %s", model))
  }
  return(res)
}

getsigwindows_quickrun=function(file = NULL, mat = NULL, threshold = 0.05, printFullOut=0,tol=10^-5,
                       diff=0,model = "MR")
{
  
  time.start <- Sys.time()
  library(MASS)
  library(hmmcov)
  
  options(scipen=999)
  
  
  #quick run
  #adjust enrichment input
  if(eadjust==T){
    #TODO: update to version that does not specify anything for enrichment, true "no distribution" scenario for enrichment
    formula = formulaB = exp_count~input_count-1
    formulaE = formulaZ = exp_count~1
  }else{
    formulaE = exp_count~input_count-1
    formula = exp_count~input_count-1
    formulaZ = exp_count~1    
  }
  
  
  #get files
  #TODO:  strengthen check for file format, maybe make sep function for it
  if(is.null(mat) & !is.null(file)){ 
    files = unlist(strsplit(files), split = ";")
  }else if(is.null(file) & !is.null(mat)){
    files = "single"
    data = mat
    if(sum(colnames(data) == "exp_count") < 1){
      stop("column names of data matrix much have one exp_count column  for experimental data")
    }else if(sum(colnames(data) == "input_count") < 1){
      print("WARNING: no input_count column found in data matrix, defaulting all formulas to exp_count~1")
      formulaE = exp_count~1
      formulaB = exp_count~1
      formulaZ = exp_count~1  
    }
  }
  
  while(fnum <= length(files)){
      # read in data
      if(length(files) > 1) data=read.table(files[fnum], header=TRUE)
      
      #TODO: add in another check for the data, maybe a check function same as above
      # set chr
      chrm = data$chromosome[1]
      
      # build covariate matrices
      mf <- model.frame(formula=formula, data=data)
      mfZ <- model.frame(formula=formulaZ, data=data)
      XB <- model.matrix(attr(mf, "terms"), data=mf)
      XZ <- model.matrix(attr(mfZ, "terms"), data=mfZ)
      
      #TODO:  this and needs to be updated to handle for no dist in E component situation
      if(eadjust==T){
        mfE <- model.frame(formula=formulaE, data=data)
        XE <- model.matrix(attr(mfE, "terms"), data=mfE)
      }else{
        XE = NULL
      }
      
      #get response data
      y <- model.response(mf)
      
      #inialization of model, choosing prop1
      if(fnum == 1){
        bgprops = c(.99, .95, .9, .8, .7)
        ll = rep(length(bgprops))
        parms = list()
        for(i in 1:length(bgprops)){
          res = run(y=y, X=XB, XE = XE, XZ = XZ, prop1 = bgprops[i],glmtype="NB", maxitEM = 5, model = model)
          ll[i] = res$ll
          parms[[i]] = res$coef
        }
        prop = bgprops[which.max(ll)]
        #recreta pp
      }
      
      #TODO:  run model, use above to initialize
      #res = run(y=y, X=XB, XE = XE, XZ = XZ, prop1 = prop, glmtype="NB", model = model, Pi, forwardback)
      res = run(y=y, X=XB, XE = XE, XZ = XZ, prop1 = .5 ,glmtype="NB", model = model)
      
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
