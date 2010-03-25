zinbagui=function(){

require(traitr)
require(gWidgets)
require(zinba)
options(guiToolkit="RGtk2")

dlg <- aDialog(items=list(
                 buildwin=choiceItem(as.integer(1), values=c(as.integer(1), as.integer(0)), tooltip="generate sample and covariate information?", editor_type="gcombobox"),
		 refinepeaks=choiceItem(as.integer(1), values=c(as.integer(1), as.integer(0)), tooltip="get refined peak boundaries?", editor_type="gcombobox"),
		 printFullOut=choiceItem(as.integer(1), values=c(as.integer(1), as.integer(0)),tooltip="print out full results in getsigwindows?",editor_type="gcombobox"),
	         method=choiceItem("mixture", values=c("mixture", "pscl"), tooltip="analysis method, mixture=mixture regression model",editor_type="gcombobox"),
		 outfile=stringItem("project name", tooltip="analysis method, mixture=mixture regression model"),	
 		 numProc=integerItem(1, tooltip="Window Size for genomic windows"),	 
		 twoBit=fileItem("", attr=list(
                                     filter=list(
				       "2bit" = list(patterns=c("*.2bit")),
                                       "All files" = list(patterns=c("*"))
                                       ))),
                 seq=fileItem("", attr=list(
                                     filter=list("Bed"=list(
                                                   patterns=c("*.bed")
                                                   ),
				       "tagAlign" = list(patterns=c("*.taf","*.tagAlign")),
                                       "All files" = list(patterns=c("*"))
                                       ))),
                 input=fileItem("none", attr=list(
                                     filter=list("Bed"=list(
                                                   patterns=c("*.bed")
                                                   ),
				       "tagAlign" = list(patterns=c("*.taf","*.tagAlign")),
                                       "All files" = list(patterns=c("*"))
                                       ))),
		 filetype=choiceItem("bed", values=c("bed","tagAlign", "bowtie"), editor_type="gcombobox"),
		 align=stringItem("Path To Align Directory", tooltip="Path to Align Directory"),	
                 filelist=stringItem("Path To filelist", tooltip="Path to filelist"),
		 winSize=integerItem(500, tooltip="Window Size for genomic windows"),
                 offset=integerItem(125, tooltip="BP offset to shift windows, window size must be a multiple of this number"),
                 extension=integerItem(200, tooltip="average length in sample library used for sequencing"),
		 formula=stringItem("", tooltip="formula for modeling background"),
		 formulaE=stringItem("exp_count~1", tooltip="formula for modeling enriched region"),
	         basecountfile=fileItem("", attr=list(
                                     filter=list(
				       "basecount" = list(patterns=c("*.basecount")),
                                       "All files" = list(patterns=c("*"))
                                       ))),
		 pWinSize=integerItem(200, tooltip='sliding window size to detect local maximums'),
		 pquant=numericItem(1, tooltip='quantile threshold for basecounts'),
	  	 threshold=numericItem(.01, tooltip='pscl method FDR threshold, lower is more stringent'),
	         peakconfidence=numericItem(.8, tooltip='mixture regression method posterior probability threshold, higher is more stringent'),
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
                         aNotebookPage(aFrame(aFrame("buildwin","refinepeaks","printFullOut","method","outfile","numProc","twoBit",label='General Options'), 

                                        aFrame("seq","input","filetype","align","filelist","winSize","offset", "extension",label='Build Window Data'),label='', horizontal=TRUE),

					aFrame(aFrame("formula", "formulaE","threshold", "peakconfidence",label='Window Analysis Options'), 

					aFrame("basecountfile","pWinSize","pquant",label='Peak Refinement Options'), label='', horizontal=TRUE), 
					label="Run Zinba"),
                         aNotebookPage(aFrame("cnvWinSize", "cnvOffset","tol", "initmethod",label=''),
                                       
                                       label="Optional Parameters")
                         )

dlg$OK_handler <- function(.) {
var=.$to_R()

do.call("run.zinba",list(filelist=var$filelist,formula=as.formula(var$formula),formulaE=as.formula(var$formulaE), outfile=var$outfile, seq=var$seq,align=var$align,input=var$input,twoBit=var$twoBit,winSize=var$winSize,offset=var$offset,cnvWinSize=var$cnvWinSize,cnvOffset=var$cnvOffset, basecountfile=var$basecountfile, threshold=var$threshold,peakconfidence=var$peakconfidence,tol=var$tol,numProc=var$numProc, buildwin=var$buildwin, pWinSize=var$pWinSize,pquant=var$pquant,refinepeaks=var$refinepeaks, printFullOut=var$printFullOut,method=var$method,initmethod=var$initmethod, diff=0, filetype=var$filetype,extension=var$extension ))
#run.zinba(var$filelist,as.formula(var$formula),as.formula(var$formulaE), var$outfile, var$seq,var$align,var$input,var$twoBit,var$winSize,var$offset,var$cnvWinSize,var$cnvOffset, var$basecountfile, var$threshold,var$peakconfidence,var$tol,var$numProc, var$buildwin, var$pWinSize,var$pquant,var$refinepeaks, var$printFullOut,var$method,var$initmethod, 0, var$filetype,var$extension )
 }
dlg$make_gui(gui_layout=view)
}
