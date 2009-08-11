#!/usr/bin/perl
use strict;
use Getopt::Long;
#use Parallel::ForkManager;

#my $TWOBIT_DIR = $ENV{TWO_BIT};

my $usage = <<'USAGE';

################ CNV-seq ################

	usage: buildPeakWindows.pl [options]

		--seq data_seq_hits.txt
		--align alignability_hits.txt
		--gdna gdna_seq_hits.txt
		--cnvarray cnvArray.tsv
		
		--window-size (default = 500 bp)
		--offset-size (default = 0)
		
		--log2-trans (default FALSE)
		--crt-trans (default FALSE)
		--db (default = hg18)

		--help

#########################################

USAGE

my ($seq_hits,$rand_hits,$gdna_hits,$cOut,%files);
my $window_size = 500;
my $offset = 0;
my $cnv_array = "none";
my $gDnaAlignTrans = "none";
my $log2Trans = "FALSE";
my $crtTrans = "FALSE";

my $result = GetOptions(
	"seq=s" => \$seq_hits,
	'window-size=i' => \$window_size,
	"offset-size=i" => \$offset,
	"log2-trans"  => sub{$log2Trans='TRUE'},
	"crt-trans"  => sub{$crtTrans='TRUE'},
	"align=s" => \$rand_hits,
	"gdna=s" => \$gdna_hits,
	"cnvarray=s" => \$cnv_array,
	"help|?" => sub{print $usage; exit}
);

die $usage unless($seq_hits);

my $nThresh = 0.5;

my @offsets = 0;
my $numOffsets = 1;
$numOffsets = int($window_size/$offset) if $offset > 0;
print STDERR "Windows will be stepped at:\n\t0bp\n";
for (my $o = 1; $o < $numOffsets; $o++){
	print STDERR "\t" . ($offset*$o) . "bp\n";
	push(@offsets,($offset*$o));
}

my (%cnvProbe,$sortCnvStarts);
if($cnv_array ne "none"){
	print STDERR "\n\nLoading cnv array data $cnv_array .....";
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
	print STDERR "FINISHED\n";
}

print STDERR "\nBuilding windows......\nUsing $seq_hits as experimental sample\n";
open(SEQ,$seq_hits) or die;
print STDERR "Counting seq hits, completed: ";
my (%count,%chrom,%maxPos);
my $pat = qr/^(\w+)\t(\d+)$/;
my $n = 0;
while(<SEQ> =~ m/$pat/){
	$n++;
	print STDERR "$n " if $n%1000000 == 0;
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

open(RAND,$rand_hits) or die;
print STDERR "FINISHED\n\nUsing $rand_hits as alignability data\nCounting align hits, completed: ";
$n = 0;
while(<RAND> =~ m/$pat/){
	$n++;
	print STDERR "$n " if $n%1000000 == 0;
	for (my $o = 0; $o <= $#offsets; $o++){
		my $window_pos = int(($2-$offsets[$o])/$window_size);
		$count{$1}->[1][$o][$window_pos]++;
	}
}close RAND;
	
open(CNV,$gdna_hits) or die;
print STDERR "FINISHED\n\nUsing $gdna_hits as input control\nCounting input hits, completed: ";
$n = 0;
while(<CNV> =~ m/$pat/){
	$n++;
	print STDERR "$n " if $n%1000000 == 0;
	for (my $o = 0; $o <= $#offsets; $o++){
		my $window_pos = int(($2-$offsets[$o])/$window_size);
		$count{$1}->[2][$o][$window_pos]++;
	}
}close CNV;

print STDERR "FINISHED\n";
my $out = $seq_hits;
$out =~ s/\..*/\_win$window_size\_/g;

#$pm->set_max_procs(25);
foreach my $chr(sort{$a<=>$b} keys %chrom){
	my $chrm = "chr" . $chr;
	for (my $o = 0; $o <= $#offsets; $o++){
		$cOut = $out . "offset" . $offsets[$o] . "bp_" . $chrm . ".temp";
		my $gcSeq = $out . "offset" . $offsets[$o] . "bp_" . $chrm . ".gcseq";
		print STDERR "\nCreating temp windows for $chrm .....";
		$files{$chr}{$offsets[$o]}{temp} = $cOut;
		$files{$chr}{$offsets[$o]}{gcSeq} = $gcSeq;
#		my $pid = $pm->start and next;
		&process_chrm($chrm,$offsets[$o],$cOut,$count{$chr},\%{$cnvProbe{$chrm}},$sortCnvStarts->{$chrm},$o) if $cnv_array ne "none";
		&process_chrm($chrm,$offsets[$o],$cOut,$count{$chr},"none",$sortCnvStarts->{$chrm},$o) if $cnv_array eq "none";
#		$pm->finish;
	}
}
#$pm->wait_all_children;
print STDERR "\n...........................FINISHED\n";

print STDERR "\nGetting seq for windows";
my %gcSeqFiles;
foreach my $chr(sort{$a<=>$b} keys %chrom){
	my $chrm = "chr" . $chr;
	for (my $o = 0; $o <= $#offsets; $o++){
		my $gcSeq = $files{$chr}{$offsets[$o]}{gcSeq};
		print STDERR "\n\tGetting gc percent for " . $files{$chr}{$offsets[$o]}{temp} ."....";
		my $start = $offsets[$o];
		$start-- if $start > 0;
		#my $end = ($window_size*$maxPos{$chr}{$offsets[$o]})+$offsets[$o]+$window_size;
		my $winOut = $files{$chr}{$offsets[$o]}{temp};
		$winOut =~ s/\..*/\.txt/g;
		my $cFileLen = $chrm . "_" . $offset . ".txt";
		`twoBitInfo /gbdb/hg18/hg18.2bit:$chrm $cFileLen`;
		open(CLEN,$cFileLen);
		my $line = <CLEN>;
		chomp($line);
		my ($cInfo,$cLength) = split(/\t/, $line);
		close CLEN; unlink($cFileLen);
		`twoBitToFa -noMask -seq=$chrm -start=$start -end=$cLength $twoBitFile $gcSeq 2> /dev/null`;
#		my $pid = $pm->start and next;
		&get_gcPerc($files{$chr}{$offsets[$o]}{gcSeq},$files{$chr}{$offsets[$o]}{temp},$winOut,$window_size);
#		$pm->finish;
	}
}
#$pm->wait_all_children;
print STDERR "\n...........................FINISHED\n";

foreach my $chr(sort{$a<=>$b} keys %chrom){
	unlink($files{$chr}{temp});
}

print STDERR "\nFINISHED\n\n";

sub process_chrm{
	my ($chrm, $offset,$cOut,$count_ref,$cnv_href,$cnvStarts,$o) = @_;
	my $cnvIndex = 0;
	my @hits = @{$count_ref->[0][$o]};
	my $numWins = @hits;
	open(OUT, ">$cOut");

	print OUT "chromosome\tstart\tend\texp_count";
	print OUT "\talign" if (defined($rand_hits) && length($rand_hits));
	print OUT "\tinput_count" if (defined($gdna_hits) && length($gdna_hits));
	print OUT "\tgdna_align_log2" if $log2Trans eq "TRUE";
	print OUT "\tgdna_align_crt" if $crtTrans eq "TRUE";
	print OUT "\tcube_cnvArray" if $cnv_array ne "none";
	print OUT "\n";

	for my $id (0..$#hits){
		
		my $start = ($window_size*$id)+$offset+1;
		my $end = $start+$window_size-1;

		my $hits_raw = $hits[$id]+0;
		print OUT "$chrm\t$start\t$end\t$hits_raw";
		
		my ($align_raw,$input_raw);
		if (defined($rand_hits) && length($rand_hits)){
			$align_raw = ${$count_ref->[1][$o]}[$id]+0;
			print OUT "\t$align_raw";
		}
		
		if (defined($gdna_hits) && length($gdna_hits)){
			$input_raw = ${$count_ref->[2][$o]}[$id]+0;
			print OUT "\t$input_raw";
		}
		
		if ($log2Trans eq "TRUE" || $crtTrans eq "TRUE"){
			my ($log2,$crt);
			$log2 = log(($input_raw+1)/($align_raw+1))/log(2) if $log2Trans eq "TRUE";
			$crt = (($input_raw+1)/($align_raw+1))**0.3 if $crtTrans eq "TRUE";
			print OUT "\t$log2" if $log2Trans eq "TRUE";
			print OUT "\t$crt" if $crtTrans eq "TRUE";
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

sub get_gcPerc {
	my ($gcFile,$tempWin,$winOut,$winSize) = @_;
	open(GC, $gcFile);
	my $header = <GC>;
	
#need to deal different size windows
	my $numLines = int($winSize/50);

	open(TEMP, $tempWin);
	open(OUT,">$winOut");
	while(<TEMP>){
		chomp;
		if($_ =~ 'chromosome'){
			print OUT "$_\tgcPerc\n";
		}else{
			my ($gcCount,$nCount) = (0,0);
#while($collectSeq == 1)$seqLen < $winSize
			for(my $g =1; $g <= $numLines; $g++){
				my $seq = <GC>;
				$gcCount += ($seq =~ tr/[G|C]//);
				$nCount += ($seq =~ tr/N//);
			}

			if(($nCount/$winSize) < 0.25){
				print OUT "$_";
				my $gcP = sprintf "%.4f", $gcCount/$winSize;
				print OUT "\t$gcP\n";
			}
		}
		
	}close TEMP;close OUT;
	unlink($tempWin);
	unlink($gcFile);
}
