#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $usage = <<'USAGE';

################ CNV-seq ################

	usage: runGetPeaks.pl [options]

		--win-file input_file_list.txt
		--threshold (default = 0.01)
		--full-model (default FALSE)
		--output-dir (default CURRENT)
		--regType (medreg or zinb)
		--medreg-method (default pfn,fn,br,lasso)
		--specify-covs (default FALSE)

#########################################

USAGE

my $medRegScript = "/Volumes/RAID_1_Striped/FAIREseq_Peak/medReg.R";
my $zinbScript = "/Volumes/RAID_1_Striped/FAIREseq_Peak/zinb.R";
my ($regType,$rScript);
my $inputFileList;
my $fullModel = 'FALSE';
my $threshold = 0.01;
my $outputDir = "";
my $mrMethod = "pfn";
my $covariates = "FALSE";

my $gcPerc = "FALSE";
my $align = "FALSE";
my $gdna = "FALSE";
my $gdnaAlignLog2 = "FALSE";
my $gdnaAlignCrt = "FALSE";
my $cnvarray = "FALSE";

my $result = GetOptions(
	"win-file=s" => \$inputFileList,
	"full-model" => sub{$fullModel = 'TRUE'},
	"threshold=f" => \$threshold,
	"regType=s" => \$regType,
	"output-dir=s" => \$outputDir,
	"medreg-method=s" => \$mrMethod,
	"specify-covs=s" => \$covariates
);
#"specify-covs" => sub{$covariates = 'TRUE'}

die $usage unless($inputFileList && ($regType ne "zinb" || $regType ne "medreg"));

$rScript = $zinbScript if $regType eq "zinb";
$rScript = $medRegScript if $regType eq "medreg";

my $countFiles = 0;
my (@covNums,@covNames);

@covNames = split(/\,/, $covariates);
for (my $c = 0; $c <= $#covNames;$c++){
    $covNames[$c] = "data\$" . $covNames[$c];
}

open(LIST,$inputFileList);
while(<LIST>){
    chomp;
    my $inputFile = $_;

    $countFiles++;
    if($countFiles == 1){
	open(INPUT, $inputFile);
	my $hLine = <INPUT>;
	chomp($hLine);
	my @header = split(/\t/, $hLine);
	#if($covariates eq "TRUE"){
	    for(my $h = 4; $h <= $#header; $h++){

		foreach my $cName (@covNames){
		    my $tempCovName = $cName;
		    $tempCovName =~ s/data\$//g;
		    if($header[$h] eq $tempCovName){
			push(@covNums, $h);
		    }
		}
	    }
	#        print STDERR "Use $header[$h] as a covariate: [1-YES|0-NO] ";
	#        my $ans = <STDIN>;
	#        chomp($ans);
	#        push(@covNums, ($h+1)) if $ans == 1;
	#        my $cv = "data\$" . $header[$h] if $ans == 1;
	#        push(@covNames, $cv) if $ans == 1;
	#    }
	#}else{
	#    for(my $h = 4; $h <= $#header; $h++){
	#	push(@covNums, ($h+1));
	#	my $cv = "data\$" . $header[$h];
	#        push(@covNames, $cv);
	#    }
	#}
    }

    my $output = $outputDir . $inputFile;
    my $model = "basicModel";
    my $term = "+";
    $term = "*" if $fullModel eq "TRUE";
    $model = "fullModel" if $fullModel eq "TRUE";
    $output =~ s/\..*/\_q$threshold\_$model\.peaks/g;

    my $logFile = $outputDir . $inputFile;
    $logFile =~ s/\..*/\.Rlog/g;
    my $errorLog = $outputDir . $inputFile;
    $errorLog =~ s/\..*/\.Rerrlog/g;

    my $scriptFile = "runMedReg_Rcode_" . $inputFile;
    open(RSCRIPT, $rScript);
    open(TEMP,">$scriptFile");

    while(<RSCRIPT>){
	if($_ =~ "__COV__"){
	    my $covs = join(",", @covNums);
	    $_ =~ s/\_\_COV\_\_/$covs/g;
	    print TEMP "$_";
	}elsif($_ =~ "__MODEL__"){
	    my $mState = join($term, @covNames);
	    $_ =~ s/\_\_MODEL\_\_/$mState/g;
	    $_ =~ s/\_\_METHOD\_\_/$mrMethod/g;
	    print TEMP "$_";
	}else{
	    print TEMP "$_";
	}
    }close RSCRIPT;
    print TEMP "\n\n";
    print TEMP "sigpeaks<-z.peaks(\"$inputFile\",threshold=$threshold);\n";
    print TEMP "write.table(sigpeaks,\"$output\",row.names=FALSE,sep=\"\\t\")\n";
    close TEMP;
    print STDERR "\nGetting peaks for $inputFile\n";
    `R --vanilla --slave < $scriptFile > $logFile 2> $errorLog`;
    unlink($scriptFile);
}close LIST;




