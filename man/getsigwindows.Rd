    \name{getsigwindows}
          \Rdversion{1.0}
          \alias{getsigwindows}
          \title{Determine Signigicantly Enriched Windows}
          \description{
            Determines the set of enriched windows using either a residual-based or mixture model
          }
          \usage{getsigwindows(file,formula,formulaE,threshold=.01,peakconfidence=.8,winout,printFullOut=0,tol=10^-5,method='pscl',initmethod, diff=0)}
          \arguments{
	    \item{file}{path to file containing window counts and covariate information (from buildwindowdata()); paths to corresponding offsets can be separated by “;” as shown in any filelist file}
	    \item{winout}{Prefix used to denote the .wins, .winlist and .peaks files that are outputted by zinba}
            \item{method}{if 'mixture' is selected, the mixture regression model approach is used (must specify formula and formulaE). If 'pscl' is specified, then only needs to specify formula.  Pscl is an adhoc method that is used on FAIRE data}
	    \item{formula}{Specifies the formula desired to model the background component of counts.  Must begin with "exp_count~"  and end with a set of desired covariates separated by "+", where the covariates could be input_count (input), gcPerc (gc percentage), exp_cnvwin_log (cnv estimate), and align_perc (alignability estimate).  Examples:  If sufficiently sequenced input is present, exp_count~input_count.  If no input is present, exp_count~gcPerc+exp_cnvwin_log+align_perc is sufficient. }
	    \item{formulaE}{Specifies the formula desired to model the enriched component.  If one is unsure about making assumptions about the relationship with covarites and enriched signal, use exp_count~1}
	    \item{peakconfidence}{posterior probability threshold  for signfiicant windows if the mixture regression model is used}
	    \item{threshold}{FDR threshold for signfiicant windows if the pscl method is used}
	    \item{tol}{Relative tolerance for convergence, default  is 10^-5}
	    \item{initmethod}{Method of intial partitioning of windows into each component prior to mixture regression model run.  usually does nto make a difference. 'count' is default, specifies taking the windows with largest counts and classifying them intially as enriched.  'quantile' runs a quantile regression onthe data using the background covariates, takes the top residuals and assigns them as enriched.  'pscl' uses the pscl method to find likely enriched regions and assigns these windows to the enriched component. } 
	    \item{printFullOut}{If set to 1, prints out the original dataset along with the posterior probabilities (default).  Otherwise, prints out only window coordinates and posterior probabilites}
}
          \seealso{
            \code{\link{save}}.
          }
          \examples{
	   #FAIRE-seq with recommended parameters, no input, multiple offsets.  
	   run.zinba(file='FAIRE/faire_chr10_win250bp_offset0bp.txt;FAIRE/faire_chr10_win250bp_offset50bp.txt',winout='FAIRE/faire',formulaE=exp_count~1,formula=exp_count~exp_cnvwin_log+gcPerc+align_perc,peakconfidence=.95,tol=10^-5,printFullOut=1,method='mixture',initmethod='count')


   
          }
          \keyword{file}
