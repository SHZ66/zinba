    \name{buildwindowdata}
          \Rdversion{1.0}
          \alias{buildwindowdata}
          \title{Generate Data Set for ZINBA given Raw Data}
          \description{
            Starting from mapped reads/raw solexa data?, reference genome in 2bit format, and appropriately formatted covariate files, generates window data for read counts and associated covariate values to use in getsigwindows().  Looks up Perl script buildPeakWindows_V3.pl in src directory of package to process raw data. Need to check to see if the command pasted is implemented properly. 
          }
          \usage{buildwindowdata(seq, align, gdna, cnvarray)}
          \arguments{
            \item{seq}{path to raw solexa/mapped reads? input file containing mapped reads}
	    \item{align}{path to file containing alignability track}
	    \item{gdna}{path to gdna file?}
            \item{cnvarray}{path to cnv file?}
          }
          \seealso{
            \code{\link{save}}.
          }
          \examples{

   
          }
          \keyword{file}
