    \name{convertmappability}
          \Rdversion{1.0}
          \alias{convertmappability}
          \title{Convert binary Mappability Files from Peakseq into Wig Format}
          \description{
            Converts the binary output from the PeakSeq mappability script to a text format that can be input into the alignAdjust function.

          }
          \usage{convertmappability(inputfile, outputfile)}
          \arguments{
	    \item{inputfile}{binary file containing a score at every base corresponding to the number of times a given k-mer occurs, originating from that base; downloaded from repository or self generated}
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
