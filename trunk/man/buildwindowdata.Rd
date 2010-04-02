    \name{buildwindowdata}
          \Rdversion{1.0}
          \alias{buildwindowdata}
          \title{Build Datasets for Mixture Regression Model}
          \description{
            Takes the raw mapped sample reads, raw mapped input reads (if available), alignability directory, and current build of the human genome and building the datasets needed to run the analysis to detect enriched regions.  
          }
          \usage{buildwindowdata(seq,input="none",align,twoBit,winSize=500,offset=0,cnvWinSize=100000,cnvOffset=0,filelist,filetype="bowtie",extension)}
          \arguments{
            \item{filelist}{a file will be created at path 'filelist' upon building of the processed data from the raw data}
	    \item{seq}{path to mapped sample reads, formatted as either 'bed', 'tagAlign', or 'bowtie'}
	    \item{align}{path to directory containing alignability files for each chromsome, obtained from alignAdjust (check how these files are generated) or downloadd from our respository.  Alignability information is specific to the uniquness threshold one used to initially filter their mapped reads and the length of the sequence tags used}
	    \item{input}{path to mapped input reads, formatted as either 'bed', 'tagAlign', or 'bowtie'.  If left blank, then defaults to 'none'}
	    \item{extension}{average length of fragments in fragment library used, typically around 200}
	    \item{twobit}{path to current build of human genome in .2bit format}
	    \item{winSize}{window size, default=500bp}
	    \item{offset}{bp offset, default is none.  If one's winSize is 500 and would like 4 equally spaced offsets, then specifying offset=125 would achieve this}
	    \item{filetype}{format of mapped sample and input reads.  'bed' is files in the standard .bed format, 'tagAlign' signifies those in .taf format, and 'bowtie' signifies mapped reads directly outputted from bowtie.  Default is 'bowtie'}
	    \item{cnvWinSize}{Size of windows used to calculate CNV activity in sample, default size is 100000}
	    \item{cnvOffset}{Offset for CNV windows, typically 2500bp is suffcient, although default is no offsets}
	  }
          \seealso{
            \code{\link{save}}.
          }
          \examples{
#Chip-Seq with recommended parameters
buildwindowdata(
   filelist='data/ctcf.list',
   seq='data/ctcfGM12878rep3chr22.taf' 
   align='align1/',
   input='data/inputGM12878rep3chr22.taf',
   twoBit='hg18.2bit',
   winSize=500,
   offset=125,
   cnvWinSize=100000,
   cnvOffset=2500,
   filetype='tagAlign',    
   extension=200)

#FAIRE-seq with recommended parameters, no input.  
buildwindowdata(
   filelist='data/faire.list',
   seq='data/faireGM12878rep1chr22.taf',
   align='align4/',
   input='none',
   twoBit='hg18.2bit',
   winSize=250,
   offset=50,
   cnvWinSize=100000,
   cnvOffset=2500,
   filetype='tagAlign',    
   extension=134)


   
          }
          \keyword{file}
