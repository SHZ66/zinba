    \name{simulatecounts}
          \Rdversion{1.0}
          \alias{buildwindowdata}
          \title{Simulate counts given covariates and fitted ZINB regression object}
          \description{
            Given a covariate matrix with p columns for predictors and a fitted ZINB regression object from zicounts or zeroinfl, the function simulates window counts.  One can also input the number of peak regions desired in the model, in addition to how large the peaks should be.  As of now, only counts can be simulated, basecounts will be the next feature.  The function outputs a two column matrix, the first column being an indicator for whether that window is a peak or not, and the second is the simulated count for that window. 
          }
          \usage{simulatecounts(a, X, numpeaks, peakcount,method)}
          \arguments{            \item{a}{fitted ZINB regression object from zicounts or zeroinfl }
	    \item{X}{nxp covariate matrix}
	    \item{numpeaks}{approximate number of peaks desired (will not be exactly this number)}
            \item{peakcount}{size of peaks}
	    \item{method}{either 'zicounts' or 'pscl'}
          }
          \seealso{
            \code{\link{save}}.
          
          \examples{

   
          }
          \keyword{file}
