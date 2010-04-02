    \name{convertmappability}
          \Rdversion{1.0}
          \alias{convertmappability}
          \title{Convert binary Mappability Files from Peakseq into Wig Format}
          \description{
            Converts the binary output from the PeakSeq mappability script into a single wig file that can be input into the alignAdjust function.  The mappability directory that you have downloaded should contain a file called "map.list", which specifies the name of each file in the directory.  To run this function properly, you need to change your working directory to the mappability directory.  For example, if your mappability directory is 'map36/', you need to run setwd('map36') before runncing convert mappability.  

          }
          \usage{convertmappability(inputfile, outputfile)}
          \arguments{
	    \item{inputfile}{Path to file in the current directory map.list, which contains the names of the chromosome specific binary files in the mappability directory. Each chromosome binary file contains a score at every base corresponding to the number of times a given k-mer occurs in the genome, originating from that base; downloaded from repository or self generated}
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
