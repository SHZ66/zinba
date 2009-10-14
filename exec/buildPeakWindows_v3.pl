#!/usr/bin/perl
use strict;
use Getopt::Long;
use FindBin qw($Bin);
BEGIN { unshift @INC,$Bin; }
#use Parallel::ForkManager;
use ForkManager_pg;
#BEGIN { push @INC,'.'; } 
#require Parallel::ForkManager;

my $usage = <<'USAGE';

################ CNV-seq ################

	usage: buildPeakWindows.pl [options]

		--seq data_seq_hits.txt
		--align alignability_dir/
		--input input_seq_hits.txt
		--rand rand_read_hits.txt
		--cnvarray cnvArray.tsv
		--cnv-expwin cnv_expSlideWin.txt
		--cnv-expcustom cnv_expCustom.txt
		--twoBit gb.2bit
		
		--window-size (default = 500 bp)
		--offset-size (default = 0)
		
		--trans-input trasnform input counts, cube root (default FALSE)
		--input_rand_log2 calc log2 ratio of input counts over rand (default FALSE)
		--gb genome build (example hg18)
		--perc-n-thresh percent of bases in window that can be N (default 0.1)
		--processes number of concurrent jobs (default 1)

		--help
#########################################

USAGE

my ($seq_hits,$cOut,%files,$gb);
my $window_size = 500;
my $offset = 0;
my $nThresh = 0.1;
my $paramFile = undef;
my $input_hits = undef;
my $twoBitFile = undef;
my $align_dir = undef;
my $cnv_array = undef;
my $cnv_exp_win = undef;
my $cnv_exp_custom = undef;
my $rand_hits = undef;
my $input_rand_log2 = "FALSE";
my $transInput = "FALSE";
my $nProcesses = 1;

my $result = GetOptions(
	"seq=s" => \$seq_hits,
	"align=s" => \$align_dir,
	"input=s" => \$input_hits,
	"cnvarray=s" => \$cnv_array,
	"cnv-expwin=s" => \$cnv_exp_win,
	"cnv-expcustom=s" => \$cnv_exp_custom,
	"rand=s" => \$rand_hits,
	"input_rand_log2" => sub{$input_rand_log2 = 'TRUE'},
	"twoBit=s" => \$twoBitFile,
	"param-file=s" => \$paramFile,
	'window-size=i' => \$window_size,
	"offset-size=i" => \$offset,
	"processes=i" => \$nProcesses,
	"perc-n-thresh=f" => \$nThresh,
	"gb=s"	=> \$gb,
	"trans-input"  => sub{$transInput='TRUE'},
	"help|?" => sub{print $usage; exit}
);

die $usage unless($seq_hits && $gb);
#my $pm = new Parallel::ForkManager($nProcesses);
my $pm = new ForkManager_pg($nProcesses);

$paramFile = $seq_hits . ".list" if(!defined($paramFile));
$paramFile =~ s/^.*\///g if(!defined($paramFile));
open(LIST,">>$paramFile");

my @offsets = 0;
my $numOffsets = 1;
$numOffsets = int($window_size/$offset) if $offset > 0;
for (my $o = 1; $o < $numOffsets; $o++){
	push(@offsets,($offset*$o));
}

open(SEQ,$seq_hits) or die;
my (%count,%chrom,%maxPos);
my $pat = qr/^(\w+)\t(\d+)$/;
while(<SEQ> =~ m/$pat/){
	unless($1 =~ '_'){
		$chrom{$1}=1;
		for (my $o = 0; $o <= $#offsets; $o++){
			my $window_pos = int(($2-$offsets[$o])/$window_size);
			$count{$1}->[0][$o][$window_pos]++;
		}
	}
}close SEQ;

if($input_hits){
	open(INPUT,$input_hits) or die;
	while(<INPUT> =~ m/$pat/){
		for (my $o = 0; $o <= $#offsets; $o++){
			my $window_pos = int(($2-$offsets[$o])/$window_size);
			$count{$1}->[1][$o][$window_pos]++;
		}
	}close INPUT;
}

if($rand_hits){
	open(RAND,$rand_hits) or die;
	while(<RAND> =~ m/$pat/){
		for (my $o = 0; $o <= $#offsets; $o++){
			my $window_pos = int(($2-$offsets[$o])/$window_size);
			$count{$1}->[2][$o][$window_pos]++;
		}
	}close RAND;
}

my (%cnvProbe,$sortCnvStarts);
if($cnv_array){
	open(CNV,$cnv_array) or die;
	while(<CNV>){
	    chomp;
	    unless($_ =~ '#'){
		my ($probe,$chrm,$start,$stop,$ratio) = split(/\t/, $_);
		$chrm =~ s/chr//g;
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
}else{
	foreach my $chrm (%chrom){
		$cnvProbe{$chrm} = 0;
		@{$sortCnvStarts->{$chrm}} = 0;
	}
}

my(%cnvExpCoords,$cnv_winSize,$cnv_offset);
my @cnv_offsets = 0;
if($cnv_exp_win){
	open(CNVEXP,$cnv_exp_win);
	while(<CNVEXP>){
		chomp;
		if($_ =~ /\#PARAMS/){
			my @head = split(/\t/, $_);
			($cnv_winSize,$cnv_offset) = ($head[1],$head[2]);
			my $cnvNumOffsets = int($cnv_winSize/$cnv_offset);
			for (my $o = 1; $o < $cnvNumOffsets; $o++){
				push(@cnv_offsets,($cnv_offset*$o));
			}
		}else{
			my @line = split(/\t/, $_);
			my $window_pos = int(((int($line[2]+$line[1])/2)-$line[3])/$cnv_winSize);
			$cnvExpCoords{$line[0]}->[$line[3]][$window_pos] = $line[6]; # Avg # Reads / Alignable Base
#			$cnvExpCoords{$line[0]}->[$line[3]][$window_pos] = $line[4]; # Total Number Reads in Window
		}
	}close CNVEXP;
}else{
	foreach my $chrm (%chrom){
		$cnvExpCoords{$chrm}->[0][0] = 0;
	}
}

my (%cnvCusExpCoords,$cnvExpStarts);
if($cnv_exp_custom){
	open(CNVEXP,$cnv_exp_custom);
	while(<CNVEXP>){
		chomp;
		unless($_ =~ /\#PARAMS/){
			my @line = split(/\t/, $_);
			my $chrm = $line[0];
			$chrm =~ s/chr//g;
			$cnvCusExpCoords{$chrm}{$line[1]} = {
				stop => $line[2],
				score => $line[6]
			};
			push(@{$cnvExpStarts->{$chrm}}, $line[1]);
		}
	}close CNVEXP;
	foreach my $chrm (keys %{$cnvExpStarts}){
		@{$cnvExpStarts->{$chrm}} = sort {$a <=> $b} @{$cnvExpStarts->{$chrm}};
	}
}else{
	foreach my $chrm (%chrom){
		$cnvCusExpCoords{$chrm}{$chrm} = 0;
		@{$cnvExpStarts->{$chrm}} = 0;
	}
}

my $out = $seq_hits;
$out =~ s/^.*\///g;
$out =~ s/\..*/\_win$window_size\_/g;

foreach my $chr(sort{$a<=>$b} keys %chrom){
	my $chrm = "chr" . $chr;
	for (my $o = 0; $o <= $#offsets; $o++){
		$cOut = $out . "offset" . $offsets[$o] . "bp_" . $chrm . ".temp" if ($align_dir || $twoBitFile);
		$cOut = $out . "offset" . $offsets[$o] . "bp_" . $chrm . ".txt" if (!defined($align_dir) && !defined($twoBitFile));
		print LIST "\#DATA\t$cOut\t$chrm\t$offsets[$o]\n" if (!defined($align_dir) && !defined($twoBitFile));
		my $gcSeq = $out . "offset" . $offsets[$o] . "bp_" . $chrm . ".gcseq";
		$files{$chr}{$offsets[$o]}{gcSeq} = $gcSeq;
		$files{$chr}{$offsets[$o]}{temp} = $cOut;
		my $pid = $pm->start and next;
		&process_chrm($chrm,$offsets[$o],$cOut,$count{$chr},$cnvExpCoords{$chrm},\%{$cnvCusExpCoords{$chr}},$cnvExpStarts->{$chr},\%{$cnvProbe{$chrm}},$sortCnvStarts->{$chrm},$o);
		$pm->finish;
	}
}
$pm->wait_all_children;

if($align_dir){
	foreach my $chr(sort{$a<=>$b} keys %chrom){
		my $chrm = "chr" . $chr;
		my $aScore = $align_dir . "ALIGN_" . $gb . "_" . $chrm . "_ADJUST.wig";
		for (my $o = 0; $o <= $#offsets; $o++){
			my $tFile = $files{$chr}{$offsets[$o]}{temp};
			my $winOut = $tFile;
			$winOut =~ s/\..*/\.txt/g if !defined($twoBitFile);
			print LIST "\#DATA\t$winOut\t$chrm\t$offsets[$o]\n" if !defined($twoBitFile);
			$winOut .= 2 if $twoBitFile;
			$files{$chr}{$offsets[$o]}{temp} = $winOut;
			my $pid = $pm->start and next;
			&get_align($aScore,$tFile,$winOut,$window_size,$offsets[$o]);
			$pm->finish;
		}
	}
	$pm->wait_all_children;
}

if($twoBitFile){
	foreach my $chr(sort{$a<=>$b} keys %chrom){
		my $chrm = "chr" . $chr;
		for (my $o = 0; $o <= $#offsets; $o++){
			my $gcSeq = $files{$chr}{$offsets[$o]}{gcSeq};
			my $start = $offsets[$o];
			$start-- if $start > 0;
			my $tempFile = $files{$chr}{$offsets[$o]}{temp};
			my $winOut = $tempFile;
			$winOut =~ s/\..*/\.txt/g;
			print LIST "\#DATA\t$winOut\t$chrm\t$offsets[$o]\n";

			my $cFileLen = $chrm . "_" . $offset . ".txt";
			my $twobitinfo_in = "$twoBitFile:$chrm";
			system(qq`echo 'library(zinba);\ntwobitinfo(infile="$twobitinfo_in",outfile=$cFileLen);\n' | R --vanilla --slave >> /dev/null 2>> /dev/null`);			
#			`twoBitInfo /gbdb/hg18/hg18.2bit:$chrm $cFileLen`;
			open(CLEN,$cFileLen);
			my $line = <CLEN>;
			chomp($line);
			my ($cInfo,$cLength) = split(/\t/, $line);
			close CLEN; unlink($cFileLen);
			system(qq`echo 'library(zinba);\ntwobittofa(chrm="$chrm",start=$start,end=$cLength,twoBitFile="$twoBitFile",gcSeq="$gcSeq");\n' | R --vanilla --slave >> /dev/null 2>> /dev/null`);
#			`twoBitToFa -noMask -seq=$chrm -start=$start -end=$cLength $twoBitFile $gcSeq 2> /dev/null`;
			my $pid = $pm->start and next;
			&get_gcPerc($gcSeq,$tempFile,$winOut,$window_size);
			$pm->finish;
		}
	}
	$pm->wait_all_children;
}close LIST;

sub process_chrm{
	my ($chrm, $offset,$cOut,$count_ref,$cnvExp_href,$cnvCusExp_href,$cnvExpStartsC,$cnv_href,$cnvStarts,$o) = @_;
	my $cnvIndex = 0;
	my @hits = @{$count_ref->[0][$o]};
	my $numWins = @hits;
	open(OUT, ">$cOut");

	print OUT "chromosome\tstart\tend\texp_count";
	print OUT "\tinput_count" if $input_hits;
	print OUT "\trand_count" if $rand_hits;
	print OUT "\tcrt_input" if $transInput eq "TRUE";
	print OUT "\tinput_rand_log2" if $input_rand_log2 eq "TRUE";
	print OUT "\texp_cnvwin_est\texp_cnvwin_log" if $cnv_exp_win;
	print OUT "\texp_cnvcustom_est\texp_cnvcustom_log" if $cnv_exp_custom;
	print OUT "\tcube_cnvArray" if $cnv_array;
	print OUT "\n";

	for my $id (0..$#hits){
		
		my $start = ($window_size*$id)+$offset+1;
		my $end = $start+$window_size-1;

		my $hits_raw = $hits[$id]+0;
		print OUT "$chrm\t$start\t$end\t$hits_raw";
		
		my $input_raw;
		if ($input_hits){
			$input_raw = ${$count_ref->[1][$o]}[$id]+0;
			print OUT "\t$input_raw";
		}
		
		my $rand_raw;
		if ($rand_hits){
			$rand_raw = ${$count_ref->[2][$o]}[$id]+0;
			print OUT "\t$rand_raw";
		}
		
		if ($transInput eq "TRUE"){
			my $crt = sprintf "%.4f", ($input_raw+1)**0.3;
			print OUT "\t$crt";
		}

		if ($input_rand_log2 eq "TRUE"){
			my $log2 = sprintf "%.4f", log(($input_raw+1)/($rand_raw+1))/log(2);
			print OUT "\t$log2";
		}

		if ($cnv_exp_win){
			my $zwin = int((($start+$end)/2));
			my ($sum,$count) = (0,0);
			for (my $co = 0; $co <= $#cnv_offsets; $co++){
				my $window_pos = int(($zwin-$cnv_offsets[$co])/$cnv_winSize);
				if(defined(${$cnvExp_href->[$cnv_offsets[$co]]}[$window_pos])){
					$sum += ${$cnvExp_href->[$cnv_offsets[$co]]}[$window_pos];
					$count++;
				}
			}				
			if($count > 0){
#				my $avg = sprintf "%.2f", ($sum/$count); # Avg Number of Total Reads in Large Windows
				my $avg = sprintf "%.2f", ($sum/$count)*(2*$window_size); # Est Number Reads Based on Avg # Reads/Alignable Bases
				my $log = sprintf "%.3f", log($avg+1); 
				print OUT "\t$avg\t$log";
			}else{
				print OUT "\t0\t0";
			}
		}
		
		if($cnv_exp_custom){
			my $zwin = int((($start+$end)/2));
			my $findOverlap = 0;
			my ($maxInt,$minInt) = ($#{$cnvExpStartsC},0);
			while( ($minInt <= $maxInt) && $findOverlap == 0){
				my $node = int(($maxInt + $minInt)/2);
				if($zwin < ${$cnvExpStartsC}[$node]){
					$maxInt = $node - 1;
				}elsif($zwin > $cnvCusExp_href->{${$cnvExpStartsC}[$node]}{stop}){
					$minInt = $node + 1;
				}elsif( $zwin >= ${$cnvExpStartsC}[$node] && $zwin <= $cnvCusExp_href->{${$cnvExpStartsC}[$node]}{stop}){
					my $est_score = sprintf "%.2f", $cnvCusExp_href->{${$cnvExpStartsC}[$node]}{score} * (2*$window_size);
					my $log = sprintf "%.3f", log($est_score+1);
					print OUT "\t$est_score\t$log";
					$findOverlap = 1;
				}
			}
			print OUT "\t0\t0" if($findOverlap == 0);
		}

		if ($cnv_array){
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
	my $currStart = 1;

	open(TEMP, $tempWin);
	open(OUT,">$winOut");
	
	while(<TEMP>){
		chomp;
		if($_ =~ 'chromosome'){
			print OUT "$_\talign_count\talign_perc\n";
		}else{
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
						$alignCount += $aScore;
					}
					my $aPerc = sprintf "%.4f", $alignCount/($winSize*2);
					print OUT "$pLine\t$alignCount\t$aPerc\n";
				}
			}else{
				print OUT "$pLine\tNA\tNA\n";
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
