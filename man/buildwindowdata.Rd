    \name{buildwindowdata}
          \Rdversion{1.0}
          \alias{buildwindowdata}
          \title{Generate Data Set for ZINBA given Raw Data}
          \description{
            Tabulates mapped reads and covariate values for non-overlapping windows of user specified length accross the genome.  Currently only can tabulate a fixed number of covariates, may be modified later to be more general.  Output from this function then can be passed to getsigwindows().  Looks up Perl script buildPeakWindows_V3.pl in exec directory of package to process raw data. Has the option of being parallelized. 
          }
          \usage{buildwindowdata(paramFile=NULL,seq=NULL, align=NULL, input=NULL, cnvarray=NULL,cnvexpwin=NULL,cnvexpcustom=NULL,rand=NULL,twoBit=NULL,winSize=500,offset=0,inputTrans=NULL,inputRandLog=NULL,percN=0.1,gb="hg18",nProcess=1)}
          \arguments{
	    \item{paramfile}{}
            \item{seq}{path to mapped sample reads}
	    \item{align}{path to directory containing alignability files for each chromsomes (check how these files are generated)}
	    \item{input}{path to input control}
	    \item{rand}{path to random reads file (simulated reads)}
            \item{cnvarray}{path to cnv file in .tsv format}
	    \item{cnvexpwin}{path to sliding window version of cnv covariate, needs to be log transformed later}
	    \item{cnvexpcustom}{path to segseq version of cnv covariate, needs to be log transformed later}
	    \item{twobit}{path to current build of human genome in .2bit format}
	    \item{winsize}{window size, default=500pb}
	    \item{offset}{bp offset}
	    \item{inputTrans}{transformation of input data (log, exp, etc)}
	    \item{inputRand}{transformation of rand data (log, exp, etc)}          
	    \item{percN}{Threshold of the number of allowable bases in window that can be "N", default is .1}
	    \item{nProcess}{number of concurrent jobs to run in parallel, default is 1}
	}
          \seealso{
            \code{\link{save}}.
          }
          \examples{

   
          }
          \keyword{file}
