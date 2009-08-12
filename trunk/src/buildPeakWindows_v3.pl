#!/usr/bin/perl
use strict;
use Getopt::Long;
#use Parallel::ForkManager;

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
		
		--trans-input trasnform input counts, cube root (default FALSE)
		--gb genome build (example hg18)

		--help
#########################################

USAGE

my ($seq_hits,$cOut,%files,$gb);
my $window_size = 500;
my $offset = 0;
my $align_thresh = 1;
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
	"gb=s"	=> \$gb,
	"trans-input"  => sub{$transInput='TRUE'},
	"help|?" => sub{print $usage; exit}
);

die $usage unless($seq_hits && $gb);

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
		$cOut = $out . "offset" . $offsets[$o] . "bp_" . $chrm . ".temp";
		my $gcSeq = $out . "offset" . $offsets[$o] . "bp_" . $chrm . ".gcseq";
		$files{$chr}{$offsets[$o]}{gcSeq} = $gcSeq;
		$files{$chr}{$offsets[$o]}{temp} = $cOut;
#		my $pid = $pm->start and next;
		&process_chrm($chrm,$offsets[$o],$cOut,$count{$chr},\%{$cnvProbe{$chrm}},$sortCnvStarts->{$chrm},$o) if $cnv_array ne "none";
		&process_chrm($chrm,$offsets[$o],$cOut,$count{$chr},"none","none",$o) if $cnv_array eq "none";
#		$pm->finish;
	}
}
#$pm->wait_all_children;

if($twoBitFile ne "none"){
	foreach my $chr(sort{$a<=>$b} keys %chrom){
		my $chrm = "chr" . $chr;

		for (my $o = 0; $o <= $#offsets; $o++){
			my $gcSeq = $files{$chr}{$offsets[$o]}{gcSeq};
			my $start = $offsets[$o];
			$start-- if $start > 0;
			my $tempFile = $files{$chr}{$offsets[$o]}{temp};
			my $winOut = $tempFile;
			$winOut .= "_GC.temp" if ($align_dir ne "none");
			$winOut =~ s/\..*/\.txt/g if ($align_dir eq "none");
			$files{$chr}{$offsets[$o]}{temp} = $winOut;

			my $cFileLen = $chrm . "_" . $offset . ".txt";
			`twoBitInfo /gbdb/hg18/hg18.2bit:$chrm $cFileLen`;
			open(CLEN,$cFileLen);
			my $line = <CLEN>;
			chomp($line);
			my ($cInfo,$cLength) = split(/\t/, $line);
			close CLEN; unlink($cFileLen);
			`twoBitToFa -noMask -seq=$chrm -start=$start -end=$cLength $twoBitFile $gcSeq 2> /dev/null`;
	#		my $pid = $pm->start and next;
			&get_gcPerc($gcSeq,$tempFile,$winOut,$window_size);
	#		$pm->finish;
		}
	}
	#$pm->wait_all_children;
}

if($align_dir ne "none"){
	foreach my $chr(sort{$a<=>$b} keys %chrom){
		my $chrm = "chr" . $chr;
		my $aScore = $align_dir . "ALIGN_" . $gb . "_" . $chrm . ".wig";
		for (my $o = 0; $o <= $#offsets; $o++){
			my $winOut = $files{$chr}{$offsets[$o]}{temp};
			$winOut =~ s/\..*/\.txt/g;
	#		my $pid = $pm->start and next;
			&get_align($aScore,$files{$chr}{$offsets[$o]}{temp},$winOut,$window_size,$offsets[$o]);
	#		$pm->finish;
		}
	}
	#$pm->wait_all_children;
}

sub process_chrm{
	my ($chrm, $offset,$cOut,$count_ref,$cnv_href,$cnvStarts,$o) = @_;
	my $cnvIndex = 0;
	my @hits = @{$count_ref->[0][$o]};
	my $numWins = @hits;
	open(OUT, ">$cOut");

	print OUT "chromosome\tstart\tend\texp_count";
	print OUT "\tinput_count" if ($input_hits ne "none");
	print OUT "\tcrt_input" if $transInput eq "TRUE";
	print OUT "\tcube_cnvArray" if $cnv_array ne "none";
	print OUT "\n";

	for my $id (0..$#hits){
		
		my $start = ($window_size*$id)+$offset+1;
		my $end = $start+$window_size-1;

		my $hits_raw = $hits[$id]+0;
		print OUT "$chrm\t$start\t$end\t$hits_raw";
		
		my $input_raw;
		if ($input_hits ne "none"){
			my $input_raw = ${$count_ref->[1][$o]}[$id]+0;
			print OUT "\t$input_raw";
		}
		
		if ($transInput eq "TRUE"){
			my $crt = ($input_raw+1)**0.3;
			print OUT "\t$crt";
		}

		if ($cnv_array ne "none"){
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
	for (0..$offset){
		my $temp = <ALIGN>;
	}

	open(TEMP, $tempWin);
	open(OUT,">$winOut");
	while(<TEMP>){
		chomp;
		if($_ =~ 'chromosome'){
			print OUT "$_\talign\n";
		}else{

#need to shift scores by 1/2 avg fragment size
			print OUT "$_";
			my $alignCount = 0;
			for (1..$winSize){
				my $aScore = <ALIGN>;
				chomp($aScore);
				if($aScore <= $align_thresh && $aScore > 0){
					$alignCount++;
				}
			}
			if(!eof(ALIGN)){
				my $aPerc = sprintf "%.4f", $alignCount/$winSize;
				print OUT "\t$aPerc\n";
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
	my ($gcs,$ncount);
	my $readLen = 0;
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
				$ncount = ($seq =~ tr/N//);
				#if($ncount/length($seq) <= $nThresh)
				$readLen += length($seq);
				if($readLen < $winSize){
					$seq =~ s/[^GC]//g;
					$gcs .= $seq;
				}elsif($readLen == $winSize){
					$seq =~ s/[^GC]//g;
					$gcs .= $seq;
					my $gcP = sprintf "%.4f", length($gcs)/$winSize;
					print OUT "$outLine\t$gcP\n";
					$gcs = "";
					$readLen = 0;
					$gcFlag = 1;
				}elsif($readLen > $winSize){
					my $tempSeq = substr($seq,($readLen-$winSize));
					my $seq = substr($seq,0,($readLen-$winSize-1));
					$seq =~ s/[^GC]//g;
					$gcs .= $seq;
					my $gcP = sprintf "%.4f", length($gcs)/$winSize;
					print OUT "$outLine\t$gcP\n";
					$readLen = length($tempSeq);
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
