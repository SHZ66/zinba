zinbagui=function(){

require(traitr)
require(gWidgets)
require(zinba)
options(guiToolkit="RGtk2")
dlg <- aDialog(items=list(
                 RunChoice=choiceItem("", values=c("Build Windows?","Run Window Analysis?", "Refine Merged Significant Regions?"), tooltip="The Zinba steps you want to run, either together or alone", editor_type="gcheckboxgroup"),
		 printFullOut=choiceItem(as.integer(1), values=c(as.integer(1), as.integer(0)),tooltip="If 0, just prints out window coordinates, filter status, and window enrichment probability columns for .wins files",editor_type="gcombobox"),
		 method=choiceItem("mixture", values=c("mixture", "pscl"), tooltip="analysis method to use, either mixture regression model or pscl method (FAIRE)",editor_type="gcombobox"),
		 outfile=stringItem("prefix for output files", tooltip="path and prefix to where .wins and .peaks are placed if analysis and refinement selected together, for example '/home/data/analysis_1'"),	
 		 numProc=integerItem(1, tooltip="How many processing cores to allocate, more processors means faster execution but more memory needed"),	 
		 twoBit=fileItem("", attr=list(
                                     filter=list(
				       "2bit" = list(patterns=c("*.2bit")),
                                       "All files" = list(patterns=c("*"))
                                       )), tooltip="Select .2bit file for analysis"),
                 seq=fileItem("", attr=list(
                                     filter=list("Bed"=list(
                                                   patterns=c("*.bed")
                                                   ),
				       "tagAlign" = list(patterns=c("*.taf","*.tagAlign")),
                                       "All files" = list(patterns=c("*"))
                                       )),tooltip="Select mapped sample reads file for analysis"),
                 input=fileItem("none", attr=list(
                                     filter=list("Bed"=list(
                                                   patterns=c("*.bed")
                                                   ),
				       "tagAlign" = list(patterns=c("*.taf","*.tagAlign")),
                                       "All files" = list(patterns=c("*"))
                                       )),tooltip="Select mapped input control reads file (if available) for analysis"),
		 filetype=choiceItem("bed", values=c("bed","tagAlign", "bowtie"), editor_type="gcombobox",tooltip="Select format of sample and input mapped read file"),
		 align=fileItem(attr=c(type="selectdir", text="Path To Align Directory")),	
                 createfilelist=stringItem("Path To filelist", tooltip="Path to filelist that will be created, holds locations of build datasets"),
		 winSize=integerItem(500, tooltip="Window Size for data, smaller is better for low signal to noise, high background data"),
                 offset=integerItem(125, tooltip="BP offset to shift windows, window size must be a multiple of this number.  Yields better sensitivity"),
                 extension=integerItem(200, tooltip="average fragment length in sample library used for sequencing reads, usually 150bp-200bp"),
                 filelist=fileItem("Path existing file list", attr=list(
                                     filter=list(
				       "file list" = list(patterns=c("*.list")),
                                       "All files" = list(patterns=c("*"))
                                       ))),
		 formula=stringItem("", tooltip="formula for modeling background"),
		 formulaE=stringItem("exp_count~1", tooltip="formula for modeling enriched region"),
		 threshold=numericItem(.01, tooltip='pscl method FDR threshold, lower is more stringent'),
	         peakconfidence=numericItem(.8, tooltip='mixture regression method posterior probability threshold, higher is more stringent'),
	         basecountfile=fileItem("", attr=list(
                                     filter=list(
				       "basecount" = list(patterns=c("*.basecount")),
                                       "All files" = list(patterns=c("*"))
                                       ))),
		 winlist=fileItem("Path to.winlist from previous window analysis", attr=list(
                                     filter=list(
				       "winlist" = list(patterns=c("*.winlist")),
                                       "All files" = list(patterns=c("*"))
                                       ))),
		 peaksfilename=stringItem("Path To file for refined peaks", tooltip="Path To file for refined peaks"),
		 pWinSize=integerItem(200, tooltip='sliding window size to detect local maximums, smaller values are more sensitive but may yield insignificant peaks'),
		 pquant=numericItem(1, tooltip='quantile threshold for basecounts, 1 means taking only the max of the entire merged region'),
		 cnvWinSize=integerItem(100000, tooltip="Window size for estimating cnv activity"),
                 cnvOffset=integerItem(2500, tooltip="BP offset to shift CNV windows, CNV window size must be a multiple of this number"),
		 tol=numericItem(.00001, tooltip='mixture regression convergence relative tolerance level'),
		 initmethod=choiceItem("count", values=c("count","quantile", "pscl"),editor_type="gcombobox"),             
                 
  ## our things
  assign.to=stringItem("df", label="Assign to:"),
  output=tableItem(attr=list(size=c(400,400)), show_label=FALSE),
  file.type=stringItem("")
  ),
title="Zinba GUI",
help_string="Select a file, then adjust parameters.")
view <- aNotebook(
                         aNotebookPage( aFrame(
						aFrame(aContainer("RunChoice" ,"printFullOut","method","outfile","numProc","twoBit"),label='General Options'), 

                                        	aFrame(aContainer("seq","input","filetype","align","createfilelist","winSize","offset", "extension"),label='Build Window Data',
							enabled_when=function(.) {
  								val <- .$to_R()$RunChoice
  								sum(val=="Build Windows?")==1
							} 
						)
					,label='', horizontal=TRUE),

					aFrame(
						aFrame(aContainer(
							aContext("filelist", context=dlg, 
								enabled_when=function(.) {
  								val <- .$to_R()$RunChoice
  								sum(val =="Build Windows?")==0
							})									
							,"formula", 
							aContext("formulaE", context=dlg, 
								enabled_when=function(.) {
  								val <- .$to_R()$method
  								val =="mixture"
							}),
							aContext("threshold", context=dlg, 
								enabled_when=function(.) {
  								val <- .$to_R()$method
  								val =="pscl"
							}),
							aContext("peakconfidence", context=dlg, 
								enabled_when=function(.) {
  								val <- .$to_R()$method
  								val =="mixture"
							})),
							label='Window Analysis Options',
							enabled_when=function(.) {
  								val <- .$to_R()$RunChoice
  								sum(val=="Run Window Analysis?")==1
							}), 
						aFrame(aContainer("basecountfile",
							aContext("winlist", context=dlg, 
								enabled_when=function(.) {
  								val <- .$to_R()$RunChoice
  								sum(val =="Build Windows?")==0 & sum(val=="Run Window Analysis?")==0
							}),
							aContext("peaksfilename", context=dlg, 
								enabled_when=function(.) {
  								val <- .$to_R()$RunChoice
  								sum(val =="Build Windows?")==0 & sum(val=="Run Window Analysis?")==0
							})
							,"pWinSize","pquant"),label='Peak Refinement Options', 
						enabled_when=function(.) {
  								val <- .$to_R()$RunChoice
  								sum(val=="Refine Merged Significant Regions?")==1
							}
					), label='', horizontal=TRUE), 
			label="Run Zinba"),
                         aNotebookPage(aFrame("cnvWinSize", "cnvOffset","tol", "initmethod",label=''),
                                       
                                       label="Optional Parameters")
                         )

dlg$OK_handler <- function(.) {
var=.$to_R()
#set up exectution options
if(sum(var$RunChoice=="Build Windows?")==1){
	listvar=var$createfilelist
	build=1
}else{
	listvar=var$filelist
	build=0
}
if(sum(var$RunChoice=="Run Window Analysis?")==1){
	analysis=1
}else{
	analysis=0
}
if(sum(var$RunChoice=="Refine Merged Significant Regions?")==1){
	refine=1
}else{
	refine=0
}

#run conditional execution
	if(build==1 & analysis==0 & refine==0){
		do.call("buildwindowdata",list(seq=var$seq,align=var$align,input=var$input,twoBit=var$twoBit,winSize=var$winSize,offset=var$offset,cnvWinSize=var$cnvWinSize,cnvOffset=var$cnvOffset,filelist=listvar,filetype=var$filetype,extension=var$extension))
	}else if(build==0 & analysis==0 & refine==1){
		do.call("getrefinedpeaks",list(winlist=var$winlist,basecountfile=var$basecountfile,bpout=paste(var$outfile,'.temp', sep=''),peakout=var$peaksfilename,twoBit=var$twoBit,winSize=500,pWinSize=var$pWinSize,pquant=var$pquant,printFullOut=var$printFullOut,peakconfidence=var$peakconfidence,threshold=var$threshold,method=var$method))
	   
	}else{
		do.call("run.zinba",list(filelist=listvar,formula=as.formula(var$formula),formulaE=as.formula(var$formulaE), outfile=var$outfile, seq=var$seq,align=var$align,input=var$input,twoBit=var$twoBit,winSize=var$winSize,offset=var$offset,cnvWinSize=var$cnvWinSize,cnvOffset=var$cnvOffset, basecountfile=var$basecountfile, threshold=var$threshold,peakconfidence=var$peakconfidence,tol=var$tol,numProc=var$numProc, buildwin=build, pWinSize=var$pWinSize,pquant=var$pquant,refinepeaks=refine, printFullOut=var$printFullOut,method=var$method,initmethod=var$initmethod, diff=0, filetype=var$filetype,extension=var$extension ))
	}
 }
dlg$make_gui(gui_layout=view)
}
