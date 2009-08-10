    \name{basecountimport}
          \Rdversion{1.0}
          \alias{basecountimport}
          \title{Import Basecounts Given Coordinate File}
          \description{
            Imports basecounts from significant windows using generated coordinate file from getsigwindows().
          }
          \usage{basecountimport(inputfile, coordfile, outputfile, chromosome='all')}
          \arguments{
            \item{inputfile}{path to a precomputed genomwide basecounts for sample}
	    \item{coordfile}{path to coordinate file, output from getsignwindows()}
	    \item{outputfile}{path to where file of extracted basecounts should be places}
            \item{chromosome}{which chromosomes in coordinate file to get basecounts for.  Default is all.}
          }
          \seealso{
            \code{\link{save}}.
          }
          \examples{

   
          }
          \keyword{file}

