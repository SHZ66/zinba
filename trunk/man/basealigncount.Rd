    \name{basealigncount}
          \Rdversion{1.0}
          \alias{basealigncount}
          \title{Get SBC (basecount) for Mapped Sample Reads}
          \description{
            This function creates the SBC (basecount) file required by run.zinba that is later used in the peak refinement step.  This should be generated prior to running run.zinba.   
          }
          \usage{basealigncount(inputfile,outputfile,twoBitFile,extension,filetype="bowtie")}
          \arguments{
	    \item{inputfile}{path to mapped sample reads if buildwin=1, formatted as either 'bed', 'tagAlign', or 'bowtie'}
	    \item{outputfile}{path to file where basecount file will be outputted}
	    \item{twobitfile}{path to current build of human genome in .2bit format}
	    \item{extension}{average length of fragments in fragment library used, typically around 200}
	    \item{filetype}{format of mapped sample and input reads.  'bed' is files in the standard .bed format, 'tagAlign' signifies those in .taf format, and 'bowtie' signifies mapped reads directly outputted from bowtie.  Default is 'bowtie'}
	  }
          \seealso{
            \code{\link{save}}.
          }
          \examples{
	   #Corresponding to FAIRE example in run.zinba:
	   basealigncount(extension=134,   filetype='tagAlign', inputfile='faire.taf', outputfile='faire.basecount' , twoBitFile='hg18.2bit' )
   
          }
          \keyword{file}
