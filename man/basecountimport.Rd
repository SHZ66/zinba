    \name{basecountimport}
          \Rdversion{1.0}
          \alias{basecountimport}
          \title{Import Basecounts Given Coordinate File}
          \description{
            Imports basecounts from a specified basecount file given a coordinate file outputted from getsigwindow().  Basecountimport determines the max for each set of significant windows generated for each offset, then discards windows where the max occurs at the window boundary or extra windows that share the same maximum.  Windows that have maxes less than 50 bp of another higher max are also discarded.  After a set of unique window maxes have been determined, d basecounts to the left and right of each max are imported, where d is currently set to 250bp by default (will change).}
          \usage{basecountimport(inputfile, coordfile, outputfile, chromosome='all')}
          \arguments{
            \item{inputfile}{path to a precomputed genomewide basecounts for sample}
	    \item{coordfile}{path to coordinate file, output from getsigwindows()}
	    \item{outputfile}{path to where file of extracted basecounts should be placed}
            \item{chromosome}{which chromosomes in coordinate file to get basecounts for.  Default is all.}
          }
          \seealso{
            \code{\link{save}}.
          }
          \examples{

   
          }
          \keyword{file}
