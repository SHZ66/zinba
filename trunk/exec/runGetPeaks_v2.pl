#!/usr/bin/perl

use strict;
use warnings;
BEGIN { push @INC,'.'; }
require Parallel::ForkManager;
use Getopt::Long;

my ($inputFileList,$bpCountFile);
my $threshold = 0.01;
my $covariates = "FALSE";
my $concurr_process = 1;
my $winsize = 0;
my $method = "zicounts";

my $result = GetOptions(
	"win-file=s" => \$inputFileList,
	"threshold=f" => \$threshold,
	"basecount_file=s" => \$bpCountFile,
	"win-size=i" => \$winsize,
	"method=s" => \$method,
	"covs=s" => \$covariates,
	"processes=i" => \$concurr_process
);

my $pm = new Parallel::ForkManager($concurr_process);

open(LIST,$inputFileList);
while(<LIST>){
    chomp;
    my ($inputFile,$chrm) = split(/\t/,$_);
    
    my $coordout = $inputFile . "_PEAK_COORDS.temp";
    my $winout = $inputFile;
    $winout =~ s/\..*/\.wins/g;
    my $peakout = $inputFile;
    $peakout =~ s/\..*/\.peaks/g;
    my $bpOut = $coordout . "_BPcount";

#	    getsigwindows(file=as.character(data_list[i,1]),covnames=covs,threshold=threshold,winout=winout,coordout=coordout,offset=(winSize/2),method=method)
#	    basecountimport(inputfile=basecountfile,coordfile=coordout,outputfile=bpOut,chromosome=data_list[i,2])
#	    peakbound(profile=bpOut,output=peakout)

    my $pid = $pm->start and next;
    &run_zinba($inputFile,$coordout,$winout,$peakout,$bpOut,$covariates,$threshold,$winsize,$method,$bpCountFile,$chrm);
    $pm->finish;
}close LIST;
$pm->wait_all_children;

sub run_zinba{
    my ($inputFile,$coordout,$winout,$peakout,$bpOut,$covs,$threshold,$winSize,$method,$bpCountFile,$chrm) = @_;
    my $off = $winSize/2;
    system(qq`echo 'library(zinba);\ngetsigwindows(file="$inputFile",covnames="$covs",threshold=$threshold,winout="$winout",coordout="$coordout",offset=$off,method=$method);\nbasecountimport(inputfile="$bpCountFile",coordfile="$coordout",outputfile="$bpOut",chromosome="$chrm");\npeakbound(profile="$bpOut",output="$peakout");\n' | R --vanilla --slave > /dev/null 2> /dev/null`);
    unlink($bpOut);
    unlink($coordout);
    return;
}


