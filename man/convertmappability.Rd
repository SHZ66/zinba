    \name{convertmappability}
          \Rdversion{1.0}
          \alias{convertmappability}
          \title{Convert binary Mappability Files from Peakseq into Wig Format}
          \description{
            Converts the binary output from the PeakSeq mappability script to a text format that can be input into the alignAdjust function.

          }
          \usage{convertmappability(inputfile, outputfile)}
          \arguments{
	    \item{inputfile}{file containing a score at every base corresponding to the number of times a given k-mer occurs, originating from that base; data is formatted as a fixedStep wiggle track ( see
UCSC Genome Browser custom tracks).  This file is outputted from convertmappability}
	    \item{outDir}{directory where results will be stored}
	    \item{twobitfile}{path to current build of human genome in .2bit format}
	    \item{inputfile}{file containg a list of paths to the b.out files from the mappabiilty scripts to be converted into text.  One path per line.}
	    \item{outputfile}{single wig file where results should be printed}
          }
	  \seealso{
            \code{\link{save}}.
          }
          \examples{
	   #Corresponding to FAIRE example in run.zinba:
	   convertmappability(inputfile='mappabilitylist.txt',  outputfile='mappability.wig')
   
          }
          \keyword{file}
