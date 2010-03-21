    \name{alignAdjust}
          \Rdversion{1.0}
          \alias{alignAdjust}
          \title{Generate Alignability Track from Peakseq Mappability}
          \description{
            Generates a score for every base pair in the genome specifying whether it is mappable. To generate this information every base is thresholded based on the degree of uniqueness required during the alignment (absolute unique would require a count of 1) and shifts it by the same amount used to shift the sequenced reads.
          }
          \usage{alignAdjust(inputfile, outDir="",twoBitFile,athresh=1,adjustsize=0)}
          \arguments{
	    \item{inputfile}{file containing a score at every base corresponding to the number of times a given k-mer occurs, originating from that base; data is formatted as a fixedStep wiggle track ( see
UCSC Genome Browser custom tracks).  This file is outputted from convertmappability}
	    \item{outDir}{directory where results will be stored}
	    \item{twobitfile}{path to current build of human genome in .2bit format}
	    \item{athresh}{uniqueness threshold, number of occurrences of a given k-mer imposed during alignment (1 means absolute uniqueness)}
	    \item{adjustsize}{adjustsize one half of the extension length used for aligned reads}
	  }
          \seealso{
            \code{\link{save}}.
          }
          \examples{
	   #Corresponding to FAIRE example in run.zinba.  For FAIRE, four or less matches to the genome allowed, so athresh=4 rather than 1
	   alignAdjust('mappability.wig', outDir="align/",twoBitFile,athresh=4,adjustsize=134/2)
   
          }
          \keyword{file}
