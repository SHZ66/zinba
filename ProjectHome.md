Welcome to the ZINBA website!  Please join our [Google Group](http://groups.google.com/group/zinba-users) to stay up to date with the latest information regarding ZINBA.

We have now  updated ZINBA to **version 2.02.03** containing a few stability improvements.

To report **new bugs/issues**, please submit an issue report on the "issues" page (google account is needed).

# Contents #
  * [Wiki](http://code.google.com/p/zinba/wiki/):
    * [Installation](http://code.google.com/p/zinba/wiki/Installation):  ZINBA installation information
    * [Running ZINBA](http://code.google.com/p/zinba/wiki/UsingZINBA): Tutorial on how to run ZINBA
    * [Post Analysis](http://code.google.com/p/zinba/wiki/PostAnalysis): What to expect after running ZINBA
  * [Downloads](http://code.google.com/p/zinba/downloads/list): Section for downloads, including the latest ZINBA package and other files
  * [Issues](http://code.google.com/p/zinba/issues/list): Report issues with ZINBA
  * [Source](http://code.google.com/p/zinba/source/browse/): Subversion access for source code

# What is ZINBA? #
ZINBA (Zero Inflated Negative Binomial Algorithm) is a computational and statistical framework used to call regions of the genome enriched for sequencing reads originating from a diverse array of biological experiments.  We collectively refer to the sequencing data derived from these experiments as DNA-seq, including FAIRE-seq, ChIP-seq, and DNAase-seq experiments.

# How does it work? #
ZINBA makes decisions on whether a region of the genome is truly enriched for sequencing reads versus background signal given a set of factors modeling each type of signal. ZINBA then looks for chromosome-wide patterns of signal and their relationships with these factors in its decision-making process.  Sets of factors can be suggested using an automated model selection process.

# How does ZINBA view sequencing data? #
ZINBA separates the genome into short windows (default 250 bp) and counts the number of reads falling into each window.  When ZINBA is making decisions on whether a particular window is enriched over background, this binned representation of the data is utilized.

# ZINBA citation #
The ZINBA manuscript has been published in [Genome Biology](http://genomebiology.com/2011/12/7/R67/abstract), in addition to a [research highlight](http://genomebiology.com/2011/12/7/120) featuring the article. Please use the following to cite ZINBA:

Rashid N, Giresi P, Sun W, Ibrahim J, Lieb J. ZINBA integrates local covariates with DNA-seq data to identify broad and narrow regions of enrichment, even within amplified genomic regions. Genome Biol 2011, 12:[r67](https://code.google.com/p/zinba/source/detail?r=67).

# Contact #
ZINBA is a collaboration between the [Department of Biostatistics](http://www.sph.unc.edu/bios/) and the [Lieb Lab](http://lieblab.bio.unc.edu/) at the University of North Carolina at Chapel Hill. For general questions and reporting problems with this website, please contact Naim Rashid (homeriq5@gmail.com) and Paul Giresi (pgiresi@gmail.com).  For bugs and other issues with ZINBA, please submit an issue on the [Issues](http://code.google.com/p/zinba/issues/list)  page.

