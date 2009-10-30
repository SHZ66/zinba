#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
BEGIN { unshift @INC,$Bin; }
use ForkManager_pg;

my ($paramFile,$bpCountFile);
my $threshold = 0.01;
my $concurr_process = 1;
my $winsize = 0;
my $method = "pscl";
my $printLog = 0;
my $pwin_size = 200;
my $peakQuant = 0.75;
my $getRefinePeaks = 1;

my $result = GetOptions(
	"param-file=s" => \$paramFile,
	"threshold=f" => \$threshold,
	"basecount_file=s" => \$bpCountFile,
	"win-size=i" => \$winsize,
	"method=s" => \$method,
	"processes=i" => \$concurr_process,
        "refine_peaks=i" => \$getRefinePeaks,
	"peakwin-size=i" => \$pwin_size,
        "peak-quant=f" => \$peakQuant, 
	"print-log=i" => \$printLog
);
my $pm = new ForkManager_pg($concurr_process);

my ($formula,$chrm,$data,$inputFile,$form,$offset,$out,$outFile);
my ($coordout,$bpout,$winout,$peakout,$stdlog,$errlog,%couts,%bpouts,$filesOffsets);
open(LIST,$paramFile);
while(<LIST>){
    chomp;
    if($_ =~ 'FORMULA'){
	($form,$formula) = split(/\t/, $_);
    }elsif($_ =~ 'OUTPUT'){
	($out,$outFile) = split(/\t/,$_);
	$coordout = $outFile;# . "_PEAK_COORDS.temp";
	$winout = $outFile . ".wins";
	$peakout = $outFile . ".peaks";
	$stdlog = $outFile . "_std.log";
	$errlog = $outFile . "_err.log";
	$bpout = $outFile; # . "_BPcount";
    }elsif($_ =~ 'DATA'){
	($data,$inputFile,$chrm,$offset) = split(/\t/,$_);
	$bpouts{$chrm} = $bpout . "_" . $chrm . "_BPcount";
	$couts{$chrm} = $coordout . "_" . $chrm . "_PEAK_COORDS.temp";
	push(@{$filesOffsets->{$chrm}}, $inputFile);
    }
}close LIST;

foreach my $chrm (keys %{$filesOffsets}){
    my $offsetFiles = join(";",@{$filesOffsets->{$chrm}});
    print STDERR "Processing $chrm\n\t", join("\n\t",@{$filesOffsets->{$chrm}}), "\n";
    my $pid = $pm->start and next;
    &run_zinba($offsetFiles,$couts{$chrm},$winout,$formula,$threshold,$winsize,$pwin_size,$peakQuant,$method,$stdlog,$errlog,$printLog,$getRefinePeaks,$bpCountFile,$bpouts{$chrm},$peakout,$chrm);
    $pm->finish;
}
$pm->wait_all_children;

#if ($getRefinePeaks == 1){
#    if ($printLog == 0){
#        system(qq`echo 'library(zinba);\nbasecountimport(inputfile="$bpCountFile",coordfile="$coordout",outputfile="$bpout");\npeakbound(bpprofile="$bpout",output="$peakout",winoffset=$win_offset);\n' | R --vanilla --slave > /dev/null 2> /dev/null`);
#    }else{
#        system(qq`echo 'library(zinba);\nbasecountimport(inputfile="$bpCountFile",coordfile="$coordout",outputfile="$bpout");\npeakbound(bpprofile="$bpout",output="$peakout",winoffset=$win_offset);\n' | R --vanilla --slave >> $stdlog 2>> $errlog`);
#    }
#    unlink($bpout);
#    unlink($coordout);
#}

sub run_zinba{
    my ($inputFile,$coordout,$winout,$formula,$threshold,$winSize,$pwinSize,$pQuant,$method,$stdLog,$errLog,$printLog,$getRefinePeaks,$bpCountFile,$bpout,$peakout,$chrm) = @_;
    if ($printLog == 0){
        system(qq`echo 'library(zinba);\ngetsigwindows(file="$inputFile",formula=$formula,threshold=$threshold,winout="$winout",coordout="$coordout",getPeakRefine=$getRefinePeaks,method="$method");\n' | R --vanilla --slave > /dev/null 2> /dev/null`);
	if ($getRefinePeaks == 1){
            system(qq`echo 'library(zinba);\nbasecountimport(inputfile="$bpCountFile",coordfile="$coordout",outputfile="$bpout",chromosome="$chrm");\npeakbound(bpprofile="$bpout",output="$peakout",winSize=$pwinSize,quantile=$pQuant);\n' | R --vanilla --slave > /dev/null 2> /dev/null`);
	    unlink($bpout);
	    unlink($coordout);
	}
    }else{
        system(qq`echo 'library(zinba);\ngetsigwindows(file="$inputFile",formula=$formula,threshold=$threshold,winout="$winout",coordout="$coordout",getPeakRefine=$getRefinePeaks,method="$method");\n' | R --vanilla --slave >> $stdLog 2>> $errLog`);
	if ($getRefinePeaks == 1){
            system(qq`echo 'library(zinba);\nbasecountimport(inputfile="$bpCountFile",coordfile="$coordout",outputfile="$bpout",chromosome="$chrm");\npeakbound(bpprofile="$bpout",output="$peakout",winSize=$pwinSize,quantile=$pQuant);\n' | R --vanilla --slave >> $stdlog 2>> $errlog`);
	    unlink($bpout);
	    unlink($coordout);
	}
    }
    return(0);
}


