    \name{getsigwindows}
          \Rdversion{1.0}
          \alias{getsigwindows}
          \title{Output Coordinate file of significant windows and data file for significant windows}
          \description{
            Takes output from buildwindowdata() and determines enriched windows, outputting a coordinate file to extract basecounts from and a file containing information about significant windows (sigpeaks matrix)
          }
          \usage{getsigwindow(file,covnames,threshold,output, offset)}
          \arguments{
            \item{file}{path to input matrix, which is output from buildsigwindows()}
	    \item{covnames}{vector containing string names of covariates from input matrix to include in model}
	    \item{threshold}{desired FDR threshold, default is .01}
            \item{output}{path to output file, where coordfile is sent}
            \item{offset}{basepair offset}
	  }
          \seealso{
            \code{\link{save}}.
          }
          \examples{
   	  getsigwindows('/zinbadocs/modelmatrix', c('gcPerc', 'align', 'cnv'), 'outputfile.txt')
          }
          \keyword{file}
