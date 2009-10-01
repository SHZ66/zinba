#!/usr/bin/perl

use strict;
use warnings;
BEGIN { push @INC,'.'; }
require Parallel::ForkManager;
use Getopt::Long;

my ($paramFile,$bpCountFile);
my $threshold = 0.01;
my $concurr_process = 1;
my $winsize = 0;
my $method = "pscl";
my $printLog = 0;

my $result = GetOptions(
	"param-file=s" => \$paramFile,
	"threshold=f" => \$threshold,
	"basecount_file=s" => \$bpCountFile,
	"win-size=i" => \$winsize,
	"method=s" => \$method,
	"processes=i" => \$concurr_process,
	"print-log=i" => \$printLog
);

my $pm = new Parallel::ForkManager($concurr_process);

my ($formula,$chrm,$data,$inputFile,$form,$offset,$out,$outFile);
my ($coordout,$winout,$peakout,$stdlog,$errlog,$bpout,$filesOffsets);
open(LIST,$paramFile);
while(<LIST>){
    chomp;
    if($_ =~ 'FORMULA'){
	($form,$formula) = split(/\t/, $_);
    }elsif($_ =~ 'OUTPUT'){
	($out,$outFile) = split(/\t/,$_);
	$coordout = $outFile . "_PEAK_COORDS.temp";
	$winout = $outFile . ".wins";
	$peakout = $outFile . ".peaks";
	$stdlog = $outFile . "_std.log";
	$errlog = $outFile . "_err.log";
	$bpout = $outFile . "_BPcount";
    }elsif($_ =~ 'DATA'){
	($data,$inputFile,$chrm,$offset) = split(/\t/,$_);
	push(@{$filesOffsets->{$chrm}}, $inputFile);
    }
}close LIST;

foreach my $chrm (keys %{$filesOffsets}){
    print STDERR "Processing $chrm\n";
    my $offsetFiles = join(";",@{$filesOffsets->{$chrm}});
    my $pid = $pm->start and next;
    &run_zinba($offsetFiles,$coordout,$winout,$peakout,$bpout,$formula,$threshold,$winsize,$method,$bpCountFile,$chrm,$stdlog,$errlog,$printLog);
    $pm->finish;
}
$pm->wait_all_children;

#\nbasecountimport(inputfile="$bpCountFile",coordfile="$coordout",outputfile="$bpOut",chromosome="$chrm");\npeakbound(profile="$bpOut",output="$peakout");

sub run_zinba{
    my ($inputFile,$coordout,$winout,$peakout,$bpOut,$formula,$threshold,$winSize,$method,$bpCountFile,$chrm,$stdLog,$errLog,$printLog) = @_;
    my $off = $winSize/2;
    if ($printLog == 0){
        system(qq`echo 'library(zinba);\ngetsigwindows(file="$inputFile",formula=$formula,threshold=$threshold,winout="$winout",coordout="$coordout",offset=$off,method="$method");\n' | R --vanilla --slave > /dev/null 2> /dev/null`);
    }else{
        system(qq`echo 'library(zinba);\ngetsigwindows(file="$inputFile",formula=$formula,threshold=$threshold,winout="$winout",coordout="$coordout",offset=$off,method="$method");\n' | R --vanilla --slave >> $stdLog 2>> $errLog`);
    }
    unlink($bpOut);
    unlink($coordout);
    return;
}


