

# Using ZINBA #
A <a href='http://www.unc.edu/~nur2/zinbaweb/test_data.tgz'>test dataset</a> has been created that corresponds to the example code given in each function used below.

# Necessary file downloads #
Before using ZINBA, you need two files before you can start your analysis, given below.

## 1. Genome Build ##
The build of the genome you are using in .2bit format, downloadable below. Also downloadable from the from the <a href=' http://hgdownload.cse.ucsc.edu/downloads.html'>UCSC genome browser</a>. NOTE: YOUR DOWNLOADED BUILD MUST BE THE SAME BUILD THAT YOUR READS WERE MAPPED TO.
|**Human**|**Mouse**|**Drosophila**|**C. Elegans**|**Yeast**|
|:--------|:--------|:-------------|:-------------|:--------|
|<a href='http://hgdownload.cse.ucsc.edu/goldenPath/hg18/bigZips/hg18.2bit'>hg18</a>|<a href='http://hgdownload.cse.ucsc.edu/goldenPath/mm8/bigZips/mm8.2bit'>mm8</a>|<a href='http://www.bios.unc.edu/~nur2/dm3.2bit'>dm3</a>| <a href='http://www.bios.unc.edu/~nur2/ce6.2bit'>ce6</a> | <a href='http://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/bigZips/sacCer2.2bit'>sacCer2</a>|
|<a href='http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit'>hg19</a>|<a href='http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit'>mm9</a>|              |<a href='http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/ce10.2bit'>ce10</a>| <a href='http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.2bit'>sacCer3</a>|


## 2. Mappability File ##
A mappability folder corresponding to your <u>sequencing read length</u> and <u>genome</u> (links to compressed files are below). These files were generated using code from <a href='www.gersteinlab.org/proj/PeakSeq'>Peakseq </a>, developed by the Gerstein Lab. This file is used to derive final mappability score in ZINBA.

|**Human**|**Mouse**|**Drosophila**|**C. Elegans**|**Yeast**|
|:--------|:--------|:-------------|:-------------|:--------|
|<a href='http://www.bios.unc.edu/~nur2/map36.tgz'>hg18 36bp</a><br><a href='http://www.bios.unc.edu/~nur2/map50.tgz'>hg18 50bp</a> <table><thead><th>         </th><th><a href='http://www.bios.unc.edu/~nur2/map36_dm3.tgz'>dm3 36bp</a><br><a href='http://www.unc.edu/~nur2/map50_dm3.tgz'>dm3 50bp</a></th><th><a href='http://www.bios.unc.edu/~nur2/map36_ce4.tgz'>ce4 36bp</a><br><a href='http://www.bios.unc.edu/~nur2/map36_ce6.tgz'>ce6 36bp</a> </th><th><a href='http://zinba.googlecode.com/files/map36_sacCer2.tgz'>sacCer2 36bp</a><br><a href='http://zinba.googlecode.com/files/map50_sacCer2.tgz'>sacCer2 50bp</a> </th></thead><tbody>
<tr><td><a href='http://www.bios.unc.edu/~nur2/map36_hg19.tgz'>hg19 36bp</a><br><a href='http://www.bios.unc.edu/~nur2/map50_hg19.tgz'>hg19 50bp</a> </td><td><a href='http://www.bios.unc.edu/~nur2/map36_mm9.tgz'>mm9 36bp</a><br><a href='http://www.bios.unc.edu/~nur2/map50_mm9.tgz'>mm9 50bp</a> </td><td>              </td><td>  <a href='http://zinba.googlecode.com/files/map36_ce10.tgz'>ce10 36bp</a><br><a href='http://zinba.googlecode.com/files/map50_ce10.tgz'>ce10 50bp</a> </td><td> <a href='http://zinba.googlecode.com/files/map36_sacCer3.tgz'>sacCer3 36bp</a><br><a href='http://zinba.googlecode.com/files/map50_sacCer3.tgz'>sacCer3 50bp</a></td></tr></tbody></table>


<h1>Preparing your files</h1>
# Unpack your downloaded mappability folder by changing to the directory it is located in and running the command below, where ## is the specific name/number corresponding to the file<br>
<pre><code>wget http://www.bios.unc.edu/~nur2/map##.tgz<br>
tar -xzvf map##.tar.gz<br>
</code></pre>
# Enter R and feed the unpacked mappability folder to the generateAlignability() function to generate your alignability directory. You do not need to compute this again unless you are running samples mapped with different alignability thresholds (athresh) or whose sequencing library has a different average fragment length (extension)<br>
<pre><code>generateAlignability(<br>
  mapdir=, #mappability directory from unpacked mappability files<br>
  outdir=, #directory for processed files, used later in analysis<br>
  athresh=, #number of hits per read allowed during mapping process<br>
  extension=, #average fragment library length<br>
  twoBitFile=, #path to downloaded genome build file in .2bit format<br>
)<br>
</code></pre>
# Run the basealigncount function to generate the basecount file needed to obtain exact peak boundaries through peak refinement. If peak refinement is not desired then this step can be skipped. Otherwise, this should be generated for each set of sample reads you are analyzing. In most cases the unrefined estimates should be sufficient.<br>
<pre><code>basealigncount(<br>
  inputfile=, #mapped sample reads<br>
  outputfile=, # output path<br>
  extension=, #average fragment library length<br>
  filetype=, #either "bed", "bowtie", or "tagAlign"<br>
  twoBitFile=, #path to downloaded genome build file in .2bit format<br>
)<br>
</code></pre>


<h1>Analyze your Data</h1>
<h2>Input Formats</h2>
Three input formats are allowed for ZINBA:<br>
<br>
<b>BED</b>: The first <a href='http://genome.ucsc.edu/FAQ/FAQformat.html#format1'>6 fields</a> of the BED format<br>
<ul><li>Chromosome (has to match naming in genome build!)<br>
</li><li>Start<br>
</li><li>Stop<br>
</li><li>Name (can be any string/character)<br>
</li><li>Score (can be set to 0 if N/A)<br>
</li><li>Strand (can be set to + if unknown)</li></ul>

<b>tagAlign</b>: Similar to the above BED format, except "Seq" replaces "Name". This format is sometimes found on the UCSC genome browser.<br>
<ul><li>Chromosome (has to match naming in genome build!)<br>
</li><li>Start<br>
</li><li>Stop<br>
</li><li>Seq (sequence of mapped read)<br>
</li><li>Score (can be set to 0)<br>
</li><li>Strand (can be set to + if unknown)</li></ul>

<b>bowtie</b>: The <a href='http://bowtie-bio.sourceforge.net/manual.shtml#default-bowtie-output'>default output</a> from the Bowtie mapping software (.hits file).  <b>ZINBA ignores the last column</b> (comma-separated list of mismatch descriptors if mismatches are detected), so in total one should have 7 columns.  Also, <b>if there are spaces in the first column of your bowtie output file</b> (corresponding to read name) this will lead to importation errors in ZINBA version <= 2.01 .<br>
<br>
In general, the chromosome name started with the prefix <b>chr</b> and the names must match how they are used in the genome build file, including case.<br>
<br>
<br>
<h2>ZINBA pipeline function</h2>
Usage for the ZINBA pipeline function<br>
<pre><code>zinba(<br>
  refinepeaks=, #refine peaks? 1 for yes, 0 for no<br>
  seq=, #path to mapped experimental reads<br>
  input=, #path to mapped input reads if available (default is "none")<br>
  filetype=, #either 'bed', 'bowtie', or 'tagAlign'<br>
  threshold=, #FDR threshold, default is 0.05<br>
  align=, #path to alignability directory<br>
  numProc=, #number of CPUs to use, must be less than max available   (default 1)<br>
  twoBit=, #path to genome build in .2bit format<br>
  outfile=, #prefix for outputted files<br>
  extension=, #average fragment library length (size selected)<br>
 <br>
  ###################<br>
  OPTIONAL PARAMETERS<br>
  ###################<br>
  <br>
  basecountfile=, #path to basecount file if refinepeaks is 1<br>
  broad=, #broad setting, TRUE or FALSE (default)<br>
  printFullOut=, #print original data with enrichment estimates, 1 for yes (more space required), 0 for no (default)<br>
  interaction=, #whether or not to considering interaction during model selection, TRUE (default) or FALSE<br>
  mode=, #either "peaks" for peak calling (default) or "CNV" for calling likely amplified CNV regions for reads in "seq" (input reads are best)<br>
  FDR= #either TRUE (default) or FALSE. If false, then uses posterior probability to threshold peaks using 1-threshold<br>
)<br>
</code></pre>

This function is optimized for speed and for ease of use, utilizing only input if it is available. A more flexible function run.zinba is available where more parameters can be specified.<br>
<br>
For ENCODE reporting purposes, the zinba() convenience function can convert posterior probability to q-values for FDR thresholding if one sets FDR=TRUE (more liberal results than the posterior probability-based results in ZINBA manuscript).<br>
<br>
<h2>FAIRE-seq example</h2>
FAIRE data example: sequencing reads are 36 bp in length, and correspond to a fragment library length of 134 base pairs. The mappability file is unpacked, unloading the map36/ folder with the mappability files contained in it. We create the alignability directory at align_athresh4_extension134/, where our reads were mapped with an alignability threshold of matching up to 4 regions in the genome (otherwise they were filtered out). Download the necessary files, unpack your mappability directory, and create a directory to place your alignability files. Open a terminal and run the following, corresponding to the test unpacked dataset (a data/ directory will be created wherever you unpacked it).<br>
<br>
<pre><code>#open a terminal and download the necessary files and test dataset.<br>
#on the mac, use curl -O instead of wget<br>
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/bigZips/hg18.2bit<br>
wget http://www.bios.unc.edu/~nur2/map36.tgz<br>
wget http://www.unc.edu/~nur2/zinbaweb/test_data.tgz<br>
<br>
#unpack your mappability folder and test datasets<br>
tar -xzvf map36.tgz<br>
tar -xzvf test_data.tgz<br>
<br>
#create your directory to hold your alignability files. You can do this anyway you please, terminal command is below<br>
mkdir align_athresh4_extension134<br>
</code></pre>

Enter R and generate your alignability directory in addition to your basecount file<br>
<pre><code>library(zinba)<br>
generateAlignability(<br>
  mapdir='map36/',<br>
  outdir='align_athresh4_extension134/',<br>
  athresh=4,<br>
  extension=134,<br>
  twoBitFile='hg18.2bit'<br>
)<br>
<br>
basealigncount(<br>
  inputfile='data/faireGM12878rep1chr22.taf',<br>
  outputfile='data/faireGM12878rep1chr22.basecount',<br>
  extension=134,<br>
  filetype='tagAlign',<br>
  twoBitFile='hg18.2bit'<br>
)<br>
<br>
</code></pre>

Now, run your zinba analysis. Here we start from our mapped sample reads and our prepared files, build the datasets needed for analysis, run the mixture regression model with the model selection specified, and then run a peak boundary refinement to capture exact peak boundaries. This function is our pipeline function, tying together the smaller functions that work on different parts of the analysis. Because input control was not available for this data, this function only considers the covariates GC content, mappability, and our estimate for local background to model the data, and utilizes the model selection procedure to select the best model for each component.<br>
<br>
<pre><code>zinba(<br>
  align='align_athresh4_extension134/',<br>
  numProc=4,<br>
  seq='data/faireGM12878rep1chr22.taf',<br>
  basecountfile='data/faireGM12878rep1chr22.basecount',<br>
  filetype="tagAlign",<br>
  outfile="data/faire",<br>
  twoBit="hg18.2bit",<br>
  extension=134,<br>
  printFullOut=1,<br>
  refinepeaks=1,<br>
  broad=F,<br>
  input="none"<br>
)<br>
<br>
<br>
</code></pre>

<h2>ChIP-seq example</h2>
ChIP-seq data example: sequencing reads are 36 bp in length, and correspond to a fragment library length of 200 base pairs. The mappability file is unpacked, unloading the map36/ folder with the mappability files contained in it. We create the alignability directory at align_athresh1_extension200/, where our reads were mapped with an alignability threshold of matching up to 1 region in the genome (otherwise they were filtered out). Download the necessary files, unpack your mappability directory, and create a directory to place your alignability files. Open a terminal and run the following, corresponding to the unpacked test dataset (a data/ directory will be created wherever you unpacked it).<br>
<br>
<pre><code>#open a terminal and download the necessary files and test dataset.<br>
#on the mac, use curl -O instead of wget<br>
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/bigZips/hg18.2bit<br>
wget http://www.bios.unc.edu/~nur2/map36.tgz<br>
wget http://www.unc.edu/~nur2/zinbaweb/test_data.tgz<br>
<br>
#unpack your mappability folder and test datasets<br>
tar -xzvf map36.tgz<br>
tar -xzvf test_data.tgz<br>
<br>
#create your directory to hold your alignability files. You can do this anyway you please, terminal command is below<br>
mkdir align_athresh1_extension200<br>
</code></pre>

Enter R and generate your alignability directory in addition to your basecount file.<br>
<br>
<pre><code>library(zinba)<br>
  generateAlignability(<br>
  mapdir='map36/',<br>
  outdir='align_athresh1_extension200/',<br>
  athresh=1,<br>
  extension=200,<br>
  twoBitFile='hg18.2bit'<br>
)<br>
<br>
basealigncount(<br>
  inputfile='data/ctcfGM12878rep3chr22.taf',<br>
  outputfile='data/ctcfGM12878rep3chr22.basecount',<br>
  extension=200,<br>
  filetype='tagAlign',<br>
  twoBitFile='hg18.2bit' )<br>
</code></pre>

Now, run your zinba analysis. Here we start from our mapped sample reads and our prepared files, build the datasets needed for analysis, run the mixture regression model considering only input control, and then run a peak boundary refinement to capture exact peak boundaries.<br>
<br>
<pre><code>zinba(<br>
  align='align_athresh1_extension200/',<br>
  numProc=4,<br>
  seq='data/ctcfGM12878rep3chr22.taf',<br>
  basecountfile='data/ctcfGM12878rep3chr22.basecount',<br>
  filetype="tagAlign",<br>
  outfile="data/ctcf",<br>
  twoBit="hg18.2bit",<br>
  extension=200,<br>
  printFullOut=1,<br>
  refinepeaks=1,<br>
  broad=F,<br>
  input='data/inputGM12878rep3chr22.taf'<br>
)<br>
</code></pre>


<h1>Custom Analysis</h1>
Please note that this section and below uses ZINBA version 2.02 or later.  For versions < 2.02, you must manually set the peakconfidence level when FDR = F, for example peakfconfidence=0.95.  While the zinba() function offers ease of use, it doesn't offer much flexibility.  With the run.zinba() function, one can specify options that will allow for<br>
<br>
<ul><li>Faster processing by selecting larger window sizes<br>
</li><li>The ability to bypass model selection if model formulas are previously known<br>
</li><li>Specifying the distance threshold for merging windows<br>
</li><li>Rerunning an analysis using an already built set of windows<br>
</li><li>Calling multiple peaks within a broad region</li></ul>

In general, the user has much more control over how their data is processed and analyzed.  With the variety of data types available, this is the best option for users who are looking for a custom analysis for their genome and data type.  We describe some of the optional parameters in more detail later on.<br>
<br>
<pre><code>run.zinba(<br>
##Main parameters  <br>
  seq= # path to mapped experimental reads       <br>
  input= # path to mapped input reads if available (default is "none")       ,<br>
  filetype= # either 'bed', 'bowtie', or 'tagAlign'  <br>
  align= # path to alignability directory       <br>
  twoBit= # path to genome build in .2bit format       <br>
  winSize= # window size, default 500 bp       <br>
  offset=  # offset distance, default 0 bp     <br>
  extension= # average fragment library length (size selected)  <br>
  basecountfile= # path to basecount file if refinepeaks is 1       <br>
  selectmodel= # Either TRUE for model selection of FALSE (default, need to specify formulas)       <br>
  threshold= # FDR threshhold, default is 0.05    <br>
  FDR= #either TRUE (default) or FALSE.   <br>
  numProc= # number of CPUs to use, must be less than max available (default 1)       <br>
  winGap= # distance threshold for significant window merging (detault 0bp),<br>
  outfile= # path prefix for outputted files      <br>
    <br>
##Optional, can ignore/not specify  <br>
  #window building  <br>
  buildwin= # 1 to build windows (default), 0 to skip  <br>
  filelist= # if buildwin=0, path to .list file for existing set of built windows (usuall in same folder)   <br>
  cnvWinSize= # CNV window size, default 100000 bp.<br>
  cnvOffset=  # CNV window offset distance, default 2500 bp      <br>
    <br>
  #If Model selection  is FALSE  <br>
  formula= # background formula         <br>
  formulaE= # enrichment formula         <br>
  formulaZ= # zero-inflation formula   <br>
    <br>
  #If Model selection is TRUE  <br>
  selectchr= # one chromosome name to use for model selection, default "chr22"<br>
  selecttype= # either the abridged version "dirty" (default) or complete version "complete"     ,<br>
  selectcovs= # vector of covariate names (characters) to consider in model selection       <br>
  interaction= # consider two or three-way covariate interactions?  Default TRUE, otherwise FALSE.  <br>
    <br>
  #peak refinement  <br>
  refinepeaks= # 1 for refinement (default), 0 otherwise  <br>
  pWinSize= # sliding window size for local maximum detection (default 200 bp)         <br>
  pquant= # read overlap quantile threshold for local maximum selection (0.0-1.0 (default)     <br>
    <br>
  #miscellaneous options  <br>
  printFullOut= # 1 for printing full form of intermediate output, 0 otherwuse  <br>
  cleanup= # TRUE to delete all intermediate files, FALSE otherwise       <br>
  tol= # mixture regression EM algorithm convergence threshold, default 10^-5      <br>
  initmethod= # initialization method, default is "count", otherwise "quantile" or "pscl"       <br>
  diff= # experimental, 1 for two-sample comparison where input = path to second sample, othewise 0 (default)<br>
)<br>
</code></pre>

Notes on each parameter:<br>
<table><thead><th> <b>Parameter</b> </th><th> Note </th></thead><tbody>
<tr><td> winSize          </td><td> Selecting a larger window size increases speed of analysis but decreases resolution and sensitivity to detect enrichment.  Consider pairing with smaller offset distances to deal with potentially bisected peak regions.  The zinba() function uses 250bp  windows by default, but 500 - 1000 bp windows may be sufficient for one's analysis. </td></tr>
<tr><td> offset           </td><td> Smaller non-zero offset distances increase sensitivity but also increase computational builden.  An offset distance of at least half the window size is usually recommended. If a window size of 500 bp is specified, an offset distance of 250 bp will result in 2 sets of windows to be built - one corresponding to 0bp offset (original windows) and another set correponding to windows shifted 250 downstream from the first set.  </td></tr>
<tr><td> cnvWinSize       </td><td> Used to estimate the local background and amplifications due to potential CNVs.  Selecting a larger CNV window size decreases the sensitivity for local background estimates to be affected by enrichment regions.  However, larger windows may be less sensitive to amplifications in local background.  Default 100kb</td></tr>
<tr><td> cnvOffset        </td><td> Smaller non-zero offset distances increase sensitivity of the local background estimate to amplifications in local background signal but may be more influenced by spikes in local enrichment</td></tr>
<tr><td> winGap           </td><td> Specifies the distance threshold in bp to merge enriched window together.  Default is 0bp, where only adjacent significant windows or those overlapping from other offsets are merged together.  The broad=TRUE option from zinba() automatically sets winGap=5000.  </td></tr>
<tr><td> FDR              </td><td>  FDR = TRUE specifies the model to use the FDR threshold rather than posterior probabilities.  This typically results in more liberal peak calls. If false, then uses posterior probability to threshold peaks using 1-threshold. </td></tr>
<tr><td> selectmodel      </td><td> Specifying select model = FALSE skips the model selection process altogether and may save a significant amount of time (if you are considering many covariates for model selection).  However, one needs to specify the exact formulas that need to be used to model each component (see formula, formulaE, and formulaZ).  These can be taken from a previous zinba run where they are printed to the screen after the model selection process has finished.  </td></tr>
<tr><td> formula          </td><td> If selectmodel = FALSE, this must be specied to model the background regions of the genome.  Covariates include GC content "gcPerc", mappability "align_perc", input control "input_count", and local background covariate "exp_cnvwin_log".  The exact form should be written as exp_count ~ cov1 + cov2 +..., for example, exp~count ~ gcPerc + input_count,  Typically one should not mix input_count and exp_cnvwin_log.  If input is avaulable, using exp_count ~ input_count results in faster processing. If you do not wish to specify covariates, you can simply use exp_count~1, which utilizes only the intercept</td></tr>
<tr><td> formulaE         </td><td> Same as the above except corresponds to enrichment component</td></tr>
<tr><td> formulaZ         </td><td> Same as the above except corresponds to zero-inflated component</td></tr>
<tr><td> selectchr        </td><td> If selectmodel=TRUE, then specify one chromsome that should be used for model selection.  Typeically this is a smaller chromsome to allow for more efficient model selection. The name of the chromsome should be exactly how it is represented in its genome build, for example chromsome 22 in hg18 is referred to as "chr22"</td></tr>
<tr><td> selecttype       </td><td> Option "dirty" fixes the formulaZ at the intercept (no covariates) and then performs model selection for the background and enrichment components.  When this process is finished, then model selection proceeds for the zero inflated component given the selected background and enriched formulas.  Used by default and is much faster than option "complete", which does model selection in all components simultaneously.</td></tr>
<tr><td> selectcovs       </td><td> Vector of covariate names to consider during model selection.  That is, all models that are tested will only consist of combinations of these covariates.  Potential covariates are "gcPerc", "align_perc", "exp_cnvwin_log", and "input_count" (see formula notes above).  If only GC content and input control are to be considered, usage would be selectcovs = c("gcPerc", "input_count"). It is not advised to mix input_count and the local background estimate covariate.  Model selection is much faster if only input control is considered (default in zinba() if input control exists).  </td></tr>
<tr><td> interaction      </td><td> If interaction=TRUE (default) then all two and three-way interaction terms between covariates specified in selectcovs are used.  If set to FALSE, then they are ignored and model selection is much faster buy may not result in the best model. </td></tr>
<tr><td> pwinSize         </td><td> Size of the sliding window to determine local maximums from read overlap data in the peak refinement steps.  Should be smaller than winSize, default is 200 bp.  Smaller sizes results in greater sensitivity to potential local maximums, but also may produce may spurious local maximums caused by small blips in enrichment.  Last window sizes smooth over smaller potential local maximums but may miss real peak regions.  </td></tr>
<tr><td> pquant           </td><td> Read overlap quantile threshold for candidate local maximums.  For example, the .9 quantile corresponds to the 90th percentile of read overlap values within a peak region, and any local maximums discovered with a height below this value are ignored.  Default value is 1.0, selecting only the global maximum in the merged region (one peak per broad region unless there is a tie </td></tr>
<tr><td> printFullOut     </td><td> Controls printing options for files from the getsignwindows function where the enrichment posterior probability is calculated.  If 1, then the original data matrix is printing out (window coordinates, window read count, covariate values) are outputted along with the estimated enrichment.  This is helpful for basic exploratory analysis and plotting.  However, this also takes up much more disk space.  Setting this parameter to 0 only prints out the window coordinates and enrichment estimates. </td></tr>
<tr><td> cleanup          </td><td>  If cleanup=TRUE, then all intermediate files found in the folder outfile_files/ is deleted after the analysis is complete, where outfile is the path specified for the outfile parameter.  This saves diskspace, but one may not getrefined peaks to try different winGaps or thresholds without running the model selection and/or mixture regression analysis again. </td></tr></tbody></table>

<h2>FAIRE-seq Example</h2>
This example mirrors that of what is used in the zinba() example and will produre the same exact output.<br>
<br>
<pre><code>run.zinba(<br>
  seq='data/faireGM12878rep1chr22.taf',  <br>
  input="none",<br>
  filetype="tagAlign",<br>
  twoBit="hg18.2bit",<br>
  winSize=250,<br>
  offset=125,<br>
  extension=134,<br>
  basecountfile='data/faireGM12878rep1chr22.basecount',<br>
  align='align_athresh4_extension134/',<br>
  selectmodel=T, <br>
  selectchr = "chr22",<br>
  selecttype = "dirty",<br>
  selectcovs = c("gcPerc", "align_perc", "exp_cnvwin_log"), <br>
  interaction= T, <br>
  threshold=0.05,<br>
  refinepeaks=1, <br>
  numProc=4,<br>
  winGap=0,   <br>
  FDR=TRUE,<br>
  outfile="data/faire",<br>
  printFullOut=1,<br>
  method="mixture"<br>
)<br>
</code></pre>

<h2>ChIP-seq Example</h2>
This example mirrors that of what is used in the zinba() example and will produre the same exact output.<br>
<br>
<pre><code>run.zinba(<br>
  seq='data/ctcfGM12878rep3chr22.taf',  <br>
  input="data/inputGM12878rep3chr22.taf",<br>
  filetype="tagAlign",<br>
  twoBit="hg18.2bit",<br>
  winSize=250,<br>
  offset=125,<br>
  extension=200,<br>
  basecountfile='data/ctcfGM12878rep3chr22.basecount',<br>
  align='align_athresh1_extension200/',<br>
  selectmodel=T,<br>
  selectchr = "chr22",<br>
  selecttype = "complete",<br>
  selectcovs = c("input_count"), <br>
  interaction= T, <br>
  threshold=0.05,<br>
  refinepeaks=1, <br>
  numProc=4,<br>
  winGap=0,   <br>
  FDR=TRUE,<br>
  outfile="data/ctcf",<br>
  printFullOut=1,<br>
  method="mixture"<br>
)<br>
</code></pre>

<h1>Decreasing ZINBA Run Time</h1>
ZINBA was designed with high performance computing clusters in mind.  That is, the ability to utilize multiple computing cores for parallel processing.  On individual desktops however, the regular ZINBA parameters may pose a burden on the user in terms of computational time and resources.  This is worsened for desktops with limited memory (less than 3 GB) and a limited number of processors (less than quad-core CPUs).  The examples below show how to run zinba on the above data quickly.  However, for the greatest peak sensitivity and model stability it is recommended to use the original parameters or the zinba() function and at least running the model selection process. Increasing the window size to 1 KB and offset to 500 bp may be appropriate for broader data and will result in even faster processing.<br>
<br>
FAIRE, using model selection but ignoring interaction ( 10x faster than default)<br>
<pre><code>run.zinba(<br>
  seq='data/faireGM12878rep1chr22.taf',  <br>
  input="none",<br>
  filetype="tagAlign",<br>
  twoBit="hg18.2bit",<br>
  winSize=500,<br>
  offset=250,<br>
  extension=134,<br>
  basecountfile='data/faireGM12878rep1chr22.basecount',<br>
  align='align_athresh4_extension134/',<br>
  selectmodel=T, <br>
  selectchr = "chr22",<br>
  selecttype = "dirty",<br>
  selectcovs = c("gcPerc", "align_perc", "exp_cnvwin_log"), <br>
  interaction= F, <br>
  threshold=0.05,<br>
  refinepeaks=1, <br>
  numProc=4,<br>
  winGap=0,   <br>
  FDR=TRUE,<br>
  outfile="data/faire",<br>
  printFullOut=1,<br>
  method="mixture"<br>
)<br>
</code></pre>

FAIRE, skipping model selection ( 60x faster than default):<br>
<pre><code>run.zinba(<br>
  seq='data/faireGM12878rep1chr22.taf',  <br>
  input="none",<br>
  filetype="tagAlign",<br>
  twoBit="hg18.2bit",<br>
  winSize=500,<br>
  offset=250,<br>
  extension=134,<br>
  basecountfile='data/faireGM12878rep1chr22.basecount',<br>
  align='align_athresh4_extension134/',<br>
  selectmodel=F, <br>
  formula= exp_count ~ gcPerc + align_perc + exp_cnvwin_log,<br>
  formulaE= exp_count ~ gcPerc + align_perc + exp_cnvwin_log,<br>
  formulaZ= exp_count ~ gcPerc + align_perc + exp_cnvwin_log,<br>
  threshold=0.05,<br>
  refinepeaks=1, <br>
  numProc=4,<br>
  winGap=0,   <br>
  FDR=TRUE,<br>
  outfile="data/faire",<br>
  printFullOut=1,<br>
  method="mixture"<br>
)<br>
</code></pre>

ChIP-seq, ignoring model selection (fastest)<br>
<pre><code>run.zinba(<br>
  seq='data/ctcfGM12878rep3chr22.taf',  <br>
  input="data/inputGM12878rep3chr22.taf",<br>
  filetype="tagAlign",<br>
  twoBit="hg18.2bit",<br>
  winSize=500,<br>
  offset=250,<br>
  extension=200,<br>
  basecountfile='data/ctcfGM12878rep3chr22.basecount',<br>
  align='align_athresh1_extension200/',<br>
  selectmodel=F,<br>
  formula= exp_count ~ input_count,<br>
  formulaE= exp_count ~ input_count,<br>
  formulaZ= exp_count ~ input_count, <br>
  threshold=0.05,<br>
  refinepeaks=1, <br>
  numProc=4,<br>
  winGap=0,   <br>
  FDR=TRUE,<br>
  outfile="data/ctcf",<br>
  printFullOut=1,<br>
  method="mixture"<br>
)<br>
</code></pre>