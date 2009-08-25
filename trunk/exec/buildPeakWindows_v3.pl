#!/usr/bin/perl
use strict;
use Getopt::Long;
#
use Parallel::ForkManager;

my $usage = <<'USAGE';

################ CNV-seq ################

	usage: buildPeakWindows.pl [options]

		--seq data_seq_hits.txt
		--align alignability_dir/
		--input input_seq_hits.txt
		--cnvarray cnvArray.tsv
		--twoBit gb.2bit
		
		--window-size (default = 500 bp)
		--offset-size (default = 0)
		--align-thresh (default = 1)
		--perc-n-thresh (default 0.1)
		
		--trans-input trasnform input counts, cube root (default FALSE)
		--gb genome build (example hg18)

		--help
#########################################

USAGE

my ($seq_hits,$cOut,%files,$gb);
my $window_size = 500;
my $offset = 0;
my $align_thresh = 1;
my $nThresh = 0.1;
my $input_hits = undef;
my $twoBitFile = undef;
my $align_dir = undef;
my $cnv_array = undef;
my $transInput = "FALSE";

my $result = GetOptions(
	"seq=s" => \$seq_hits,
	"align=s" => \$align_dir,
	"input=s" => \$input_hits,
	"cnvarray=s" => \$cnv_array,
	"twoBit=s" => \$twoBitFile,
	'window-size=i' => \$window_size,
	"offset-size=i" => \$offset,
	"align-thresh=i" => \$align_thresh,
	"perc-n-thresh=f" => \$nThresh,
	"gb=s"	=> \$gb,
	"trans-input"  => sub{$transInput='TRUE'},
	"help|?" => sub{print $usage; exit}
);

die $usage unless($seq_hits && $gb);

#
my $pm = new Parallel::ForkManager(8);

my $dataFile_list = $seq_hits . ".list";
open(LIST,">$dataFile_list");

my @offsets = 0;
my $numOffsets = 1;
$numOffsets = int($window_size/$offset) if $offset > 0;
for (my $o = 1; $o < $numOffsets; $o++){
	push(@offsets,($offset*$o));
}

my (%cnvProbe,$sortCnvStarts);
if(defined($cnv_array)){
	open(CNV,$cnv_array) or die;
	while(<CNV>){
	    chomp;
	    unless($_ =~ '#'){
		my ($probe,$chrm,$start,$stop,$ratio) = split(/\t/, $_);
		$cnvProbe{$chrm}{$probe} = {
		    start   =>  $start,
		    stop    =>  $stop,
		};
		$cnvProbe{$chrm}{$probe}{sum} += $ratio;
		$cnvProbe{$chrm}{$probe}{count}++;
	    }
	}close CNV;
	foreach my $chrm (keys %cnvProbe){
		@{$sortCnvStarts->{$chrm}} = sort {  $cnvProbe{$chrm}{$a}{start} <=> $cnvProbe{$chrm}{$b}{start} } keys %{$cnvProbe{$chrm}};
		for (my $pID = 0;$pID <= $#{$sortCnvStarts->{$chrm}}; $pID++){
			if($pID > 0){
				if(($cnvProbe{$chrm}{${$sortCnvStarts->{$chrm}}[$pID-1]}{stop}+$window_size) >= $cnvProbe{$chrm}{${$sortCnvStarts->{$chrm}}[$pID]}{start}){
					$cnvProbe{$chrm}{${$sortCnvStarts->{$chrm}}[$pID-1]}{sum} += $cnvProbe{$chrm}{${$sortCnvStarts->{$chrm}}[$pID]}{sum};
					$cnvProbe{$chrm}{${$sortCnvStarts->{$chrm}}[$pID-1]}{count} += $cnvProbe{$chrm}{${$sortCnvStarts->{$chrm}}[$pID]}{count};
					$cnvProbe{$chrm}{${$sortCnvStarts->{$chrm}}[$pID-1]}{stop} = $cnvProbe{$chrm}{${$sortCnvStarts->{$chrm}}[$pID]}{stop};
					delete($cnvProbe{$chrm}{${$sortCnvStarts->{$chrm}}[$pID]});
					splice(@{$sortCnvStarts->{$chrm}},$pID,1);
				}
			}
		}
	
		foreach my $pID (@{$sortCnvStarts->{$chrm}}){
			$cnvProbe{$chrm}{$pID}{avg} = $cnvProbe{$chrm}{$pID}{sum}/$cnvProbe{$chrm}{$pID}{count};
			delete($cnvProbe{$chrm}{$pID}{sum});
			delete($cnvProbe{$chrm}{$pID}{count});
		}
	}
}
open(SEQ,$seq_hits) or die;
my (%count,%chrom,%maxPos);
my $pat = qr/^(\w+)\t(\d+)$/;
my $n = 0;
while(<SEQ> =~ m/$pat/){
	$n++;
	unless($1 =~ '_'){
		$chrom{$1}=1;
		for (my $o = 0; $o <= $#offsets; $o++){
			my $window_pos = int(($2-$offsets[$o])/$window_size);
			if(!defined($maxPos{$1}{$offsets[$o]}) && !exists($maxPos{$1}{$offsets[$o]})){
				$maxPos{$1}{$offsets[$o]} = 0;
			}
			$maxPos{$1}{$offsets[$o]} = $window_pos if $maxPos{$1}{$offsets[$o]} < $window_pos;
			$count{$1}->[0][$o][$window_pos]++;
		}
	}
}close SEQ;

if(defined($input_hits)){
	open(INPUT,$input_hits) or die;
	$n = 0;
	while(<INPUT> =~ m/$pat/){
		$n++;
		print STDERR "$n " if $n%1000000 == 0;
		for (my $o = 0; $o <= $#offsets; $o++){
			my $window_pos = int(($2-$offsets[$o])/$window_size);
			$count{$1}->[1][$o][$window_pos]++;
		}
	}close INPUT;
}

my $out = $seq_hits;
$out =~ s/\..*/\_win$window_size\_/g;

#$pm->set_max_procs(25);
foreach my $chr(sort{$a<=>$b} keys %chrom){
	my $chrm = "chr" . $chr;
	for (my $o = 0; $o <= $#offsets; $o++){
		$cOut = $out . "offset" . $offsets[$o] . "bp_" . $chrm . ".temp" if (defined($align_dir) && defined($twoBitFile));
		$cOut = $out . "offset" . $offsets[$o] . "bp_" . $chrm . ".txt" if (!defined($align_dir) && !defined($twoBitFile));
		print LIST "$cOut\t$chrm\n" if (!defined($align_dir) && !defined($twoBitFile));
		my $gcSeq = $out . "offset" . $offsets[$o] . "bp_" . $chrm . ".gcseq";
		$files{$chr}{$offsets[$o]}{gcSeq} = $gcSeq;
		$files{$chr}{$offsets[$o]}{temp} = $cOut;
#
my $pid = $pm->start and next;
		&process_chrm($chrm,$offsets[$o],$cOut,$count{$chr},\%{$cnvProbe{$chrm}},$sortCnvStarts->{$chrm},$o) if defined($cnv_array);
		&process_chrm($chrm,$offsets[$o],$cOut,$count{$chr},"none","none",$o) if !defined($cnv_array);
#
$pm->finish;
	}
}
#
$pm->wait_all_children;


if(defined($align_dir)){
	foreach my $chr(sort{$a<=>$b} keys %chrom){
		my $chrm = "chr" . $chr;
		my $aScore = $align_dir . "ALIGN_" . $gb . "_" . $chrm . ".wig";
		for (my $o = 0; $o <= $#offsets; $o++){
			my $tFile = $files{$chr}{$offsets[$o]}{temp};
			my $winOut = $tFile;
			$winOut =~ s/\..*/\.txt/g if !defined($twoBitFile);
			print LIST "$winOut\t$chrm\n" if !defined($twoBitFile);
			$winOut .= 2 if defined($twoBitFile);
			$files{$chr}{$offsets[$o]}{temp} = $winOut;
#
my $pid = $pm->start and next;
			&get_align($aScore,$tFile,$winOut,$window_size,$offsets[$o]);
#
$pm->finish;
		}
	}
#
$pm->wait_all_children;
}

if(defined($twoBitFile)){
	foreach my $chr(sort{$a<=>$b} keys %chrom){
		my $chrm = "chr" . $chr;
		for (my $o = 0; $o <= $#offsets; $o++){
			my $gcSeq = $files{$chr}{$offsets[$o]}{gcSeq};
			my $start = $offsets[$o];
			$start-- if $start > 0;
			my $tempFile = $files{$chr}{$offsets[$o]}{temp};
			my $winOut = $tempFile;
			$winOut =~ s/\..*/\.txt/g;
			print LIST "$winOut\t$chrm\n";

			my $cFileLen = $chrm . "_" . $offset . ".txt";
			`twoBitInfo /gbdb/hg18/hg18.2bit:$chrm $cFileLen`;
			open(CLEN,$cFileLen);
			my $line = <CLEN>;
			chomp($line);
			my ($cInfo,$cLength) = split(/\t/, $line);
			close CLEN; unlink($cFileLen);
			`twoBitToFa -noMask -seq=$chrm -start=$start -end=$cLength $twoBitFile $gcSeq 2> /dev/null`;
#
my $pid = $pm->start and next;
			&get_gcPerc($gcSeq,$tempFile,$winOut,$window_size);
#
$pm->finish;
		}
	}
#
$pm->wait_all_children;
}
close LIST;

sub process_chrm{
	my ($chrm, $offset,$cOut,$count_ref,$cnv_href,$cnvStarts,$o) = @_;
	my $cnvIndex = 0;
	my @hits = @{$count_ref->[0][$o]};
	my $numWins = @hits;
	open(OUT, ">$cOut");

	print OUT "chromosome\tstart\tend\texp_count";
	print OUT "\tinput_count" if defined($input_hits);
	print OUT "\tcrt_input" if $transInput eq "TRUE";
	print OUT "\tcube_cnvArray" if defined($cnv_array);
	print OUT "\n";

	for my $id (0..$#hits){
		
		my $start = ($window_size*$id)+$offset+1;
		my $end = $start+$window_size-1;

		my $hits_raw = $hits[$id]+0;
		print OUT "$chrm\t$start\t$end\t$hits_raw";
		
		my $input_raw;
		if (defined($input_hits)){
			my $input_raw = ${$count_ref->[1][$o]}[$id]+0;
			print OUT "\t$input_raw";
		}
		
		if ($transInput eq "TRUE"){
			my $crt = ($input_raw+1)**0.3;
			print OUT "\t$crt";
		}

		if (defined($cnv_array)){
			my $cnvData;
			my $winPos = int((((${$cnv_href}{${$cnvStarts}[$cnvIndex]}{stop}+${$cnv_href}{${$cnvStarts}[$cnvIndex]}{start})/2)-$offset)/$window_size);
			if($id < $winPos){
				if ($cnvIndex == 0){
					$cnvData = ${$cnv_href}{${$cnvStarts}[$cnvIndex]}{avg};
				}else{
					my $prevWinPos = int((((${$cnv_href}{${$cnvStarts}[$cnvIndex-1]}{stop}+${$cnv_href}{${$cnvStarts}[$cnvIndex-1]}{start})/2)-$offset)/$window_size);
					my $diffVal = (${$cnv_href}{${$cnvStarts}[$cnvIndex]}{avg} - ${$cnv_href}{${$cnvStarts}[$cnvIndex-1]}{avg})/($winPos-$prevWinPos);
					$cnvData = ${$cnv_href}{${$cnvStarts}[$cnvIndex-1]}{avg} + (abs($diffVal)*($id-$prevWinPos)) if $diffVal >= 0;
					$cnvData = ${$cnv_href}{${$cnvStarts}[$cnvIndex-1]}{avg} - (abs($diffVal)*($id-$prevWinPos)) if $diffVal < 0;
				}
			}elsif($id == $winPos){
				$cnvData = ${$cnv_href}{${$cnvStarts}[$cnvIndex]}{avg};
			}elsif($id > $winPos){
				if ($cnvIndex < $#{$cnvStarts}){
					$cnvIndex++;
					my $winPos = int((((${$cnv_href}{${$cnvStarts}[$cnvIndex]}{stop}+${$cnv_href}{${$cnvStarts}[$cnvIndex]}{start})/2)-$offset)/$window_size);
					my $prevWinPos = int((((${$cnv_href}{${$cnvStarts}[$cnvIndex-1]}{stop}+${$cnv_href}{${$cnvStarts}[$cnvIndex-1]}{start})/2)-$offset)/$window_size);
					my $dist = int((${$cnv_href}{${$cnvStarts}[$cnvIndex]}{start} - ${$cnv_href}{${$cnvStarts}[$cnvIndex-1]}{stop})/$window_size);
					$dist = 1 if $dist < 1;
					my $diffVal = (${$cnv_href}{${$cnvStarts}[$cnvIndex]}{avg} - ${$cnv_href}{${$cnvStarts}[$cnvIndex-1]}{avg})/$dist;
					$cnvData = ${$cnv_href}{${$cnvStarts}[$cnvIndex-1]}{avg} + (abs($diffVal)*($id-$prevWinPos)) if $diffVal >= 0;
					$cnvData = ${$cnv_href}{${$cnvStarts}[$cnvIndex-1]}{avg} - (abs($diffVal)*($id-$prevWinPos)) if $diffVal < 0;
				}elsif($cnvIndex == $#{$cnvStarts}){
					$cnvData = ${$cnv_href}{${$cnvStarts}[$cnvIndex]}{avg};
				}
			}
	
			my $rawData = sprintf "%.5f", (2**$cnvData)**3;
			print OUT "\t$rawData";
		}

		print OUT "\n";

	}close OUT;
}

sub get_align {
	my ($alignFile,$tempWin,$winOut,$winSize,$offset) = @_;
	open(ALIGN,$alignFile);
	my $alignHeader = <ALIGN>;
	my ($ahFix,$ahChrm,$ahstart,$ahstep) = split(/ /, $alignHeader);
	$ahstart =~ s/start\=//g;
	my $currStart = $ahstart;

	open(TEMP, $tempWin);
	open(OUT,">$winOut");
	
	while(<TEMP>){
		chomp;
		if($_ =~ 'chromosome'){
			print OUT "$_\talign\n";
		}else{
#need to shift scores by 1/2 avg fragment size
			my $pLine = $_;
			my @line = split(/\t/,$_);
			if(!eof(ALIGN)){
				while($line[1] > $currStart){
					my $temp = <ALIGN>;
					$currStart++;
				}

				if($line[2] < $currStart){
					print OUT "$pLine\tNA\n";
				}elsif($line[1] == $currStart){
					my $alignCount = 0;
					while($currStart <= $line[2]){
						my $aScore = <ALIGN>;
						chomp($aScore);
						$currStart++;
						if($aScore <= $align_thresh && $aScore > 0 && !eof(ALIGN)){
							$alignCount++;
						}
					}
					my $aPerc = sprintf "%.4f", $alignCount/$winSize;
					print OUT "$pLine\t$aPerc\n";
				}
			}else{
				print OUT "$pLine\tNA\n";
			}
		}
		
	}close TEMP;close OUT;
	close ALIGN;
	unlink($tempWin);
}

sub get_gcPerc {
	my ($gcFile,$tempWin,$winOut,$winSize) = @_;
	open(GC,$gcFile);
	my $gcHeader = <GC>;
	open(TEMP, $tempWin);
	open(OUT,">$winOut");
	my $gcs;
	my $readLen = 0;
	my $ncount = 0;
	while(<TEMP>){
		chomp;
		if($_ =~ 'chromosome'){
			print OUT "$_\tgcPerc\n";
		}else{
			my $outLine = $_;
			my $gcFlag = 0;
			while($gcFlag == 0){
				my $seq = <GC>;
				chomp($seq);
				$readLen += length($seq);
				if($readLen == $winSize || eof(GC)){
					$ncount += ($seq =~ tr/N//);
					$seq =~ s/[^GC]//g;
					$gcs .= $seq;
					if($ncount/$winSize < $nThresh){
						my $gcP = sprintf "%.4f", length($gcs)/$winSize;
						print OUT "$outLine\t$gcP\n";
					}
					$gcs = "";
					$readLen = 0;
					$ncount = 0;
					$gcFlag = 1;
				}elsif($readLen < $winSize){
					$ncount += ($seq =~ tr/N//);
					$seq =~ s/[^GC]//g;
					$gcs .= $seq;
				}elsif($readLen > $winSize){
					my $tempSeq = substr($seq,($readLen-$winSize));
					my $seq = substr($seq,0,($readLen-$winSize-1));
					$ncount += ($seq =~ tr/N//);
					$seq =~ s/[^GC]//g;
					$gcs .= $seq;
					if($ncount/$winSize < $nThresh){
						my $gcP = sprintf "%.4f", length($gcs)/$winSize;
						print OUT "$outLine\t$gcP\n";
					}
					$readLen = length($tempSeq);
					$ncount = ($tempSeq =~ tr/N//);
					$tempSeq =~ s/[^GC]//g;
					$gcs = $tempSeq;
					$gcFlag = 1;
				}
			}
		}
	}close TEMP;close OUT;
	close GC;
	unlink($tempWin);unlink($gcFile);
}
