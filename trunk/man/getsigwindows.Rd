    \name{getsigwindows}
          \Rdversion{1.0}
          \alias{getsigwindows}
          \title{Output Coordinate file of significant windows and data file for significant windows}
          \description{
            Takes output from buildwindowdata() and determines enriched windows, outputting a coordinate file to extract basecounts from and a file containing information about significant windows (sigpeaks matrix)
          }
          \usage{getsigwindows(file,formula,threshold=.01,peakconfidence=.8,priorpeakprop=.15,winout,coordout,offset=0,tol=10^-5,getPeakRefine=1,method='pscl'){
}
          \arguments{
            \item{file}{path to input matrix, which is output from buildsigwindows()}
            \item{formula}{model formula, of form response~covariate1+covariate2...}
	    \item{threshold}{if using ZINB regression, this specifies desired FDR threshold, default is 0.01}
	    \item{peakconfidence}{if using mixture model method, this specifies desired posterior peak probability threshold, default is 0.8}
	    \item{priorpeakprop}{if using mixture model method, this specifies starting value of the proportion of enriched windows in the data, default is 0.15}
            \item{winout}{path for output file containing portion of input data set that are selected as significantly enriched}
            \item{coordout}{path for output file containing start and stop coordinates of significant windows, used by basecountimport()}	
	    \item{offset}{unused}
	    \item{tol}{relative tolerance level for convergence of mixture model, default is 10^-5}
   	    \item{getPeakRefine}{specifies whether coordinates of significant windows should be printed for peakbound(), otherwise user can use their own methods to get exact peak boundaries.  Default is TRUE (1)}
	    \item{method}{Specifies which model should be used to select enriched windows.  Method='pscl' select the ZINB regression model, and method='mixture' uses the three component mixture model}
	  }
          \seealso{
            \code{\link{save}}.
          }
          \examples{

          }
          \keyword{file}
