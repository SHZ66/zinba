#!/usr/bin/perl
use strict;
use Getopt::Long;
BEGIN { push @INC,'.'; }
#use Parallel::ForkManager;
require Parallel::ForkManager;

my $usage = <<'USAGE';

################ CNV-seq ################

	usage: buildCNVWindows.pl [options]

		--seq data_seq_hits.txt
		--align-dir align_dir/
		--window_file cnv_coord.txt
		--window-size (default = 25000 bp)
		--offset-size (default = 0)
		--gb genome build (default hg18)
		--processes num parallel processes (default 1)
		
#########################################

USAGE

my ($seq_hits,$cOut,%files);
my $window_size = 25000;
my $window_file = undef;
my $alignDir = undef;
my $offset = 0;
my $twoBitFile = undef;
my $nProcesses = 1;
my $gb = "hg18";

my $result = GetOptions(
	"seq=s" => \$seq_hits,
	"align-dir=s" => \$alignDir,
	"window_file=s" => \$window_file,
	'window-size=i' => \$window_size,
	"offset-size=i" => \$offset,
	"gb=s" => \$gb,
	"processes=i" => \$nProcesses,
);

die $usage unless($seq_hits);
my $pm = new Parallel::ForkManager($nProcesses);

my @offsets = 0;
my $numOffsets = 1;
my (%cnvCoords,$cnvStarts);
if(!defined($window_file)){
	$numOffsets = int($window_size/$offset) if $offset > 0;
	for (my $o = 1; $o < $numOffsets; $o++){
		push(@offsets,($offset*$o));
	}
}else{
	open(CNV,$window_file);
	while(<CNV>){
		chomp;
		my ($chrm,$start,$stop) = split(/\t/, $_);
		$cnvCoords{$chrm}{$start} = {
			stop => $stop,
			count => 0
		};
		push(@{$cnvStarts->{$chrm}}, $start);
	}close CNV;
	
	foreach my $chrm (keys %{$cnvStarts}){
		@{$cnvStarts->{$chrm}} = sort {$a <=> $b} @{$cnvStarts->{$chrm}};
	}
}

open(SEQ,$seq_hits) or die;
my (%count,%chrom,$reads);
my $pat = qr/^(\w+)\t(\d+)$/;
while(<SEQ> =~ m/$pat/){
	unless($1 =~ '_'){
		$chrom{$1}=1;
		if(!defined($window_file)){
			for (my $o = 0; $o <= $#offsets; $o++){
				my $window_pos = int(($2-$offsets[$o])/$window_size);
				$count{$1}->[0][$o][$window_pos]++;
			}
		}else{
			push(@{$reads->{$1}},$2);
		}
	}
}close SEQ;

my $cnv_wins = $seq_hits;
$cnv_wins =~ s/\..*//g;
$cnv_wins .= $window_file if ($window_file);
$cnv_wins =~ s/\..*/\.cnvs/g if ($window_file);
$cnv_wins .= "_win" . $window_size . "_offset" . $offset . ".cnvs" if (!defined($window_file));

my $out = $seq_hits;
$out =~ s/\..*//g;

my @delFiles = ("tempHead.txt");
open(TEMP,">tempHead.txt");
print TEMP "#PARAMS\t$window_size\t$offset\n" if (!defined($window_file));
print TEMP "#PARAMS\tcustom\n" if (defined($window_file));
close TEMP;
my $catCMD = "cat tempHead.txt ";

foreach my $chr(sort{$a<=>$b} keys %chrom){
	my $chrm = "chr" . $chr;
	if(defined($window_file)){
		$cOut = $out . "_" . $chrm . ".temp";
		$files{$chrm}{0}{temp} = $cOut;
		my $pid = $pm->start and next;
		&bst_readsTOcnv(\@{$cnvStarts->{$chr}},\%{$cnvCoords{$chr}},\@{$reads->{$chr}},$chrm,$cOut);
		$pm->finish;
	}else{
		for (my $o = 0; $o <= $#offsets; $o++){
			$cOut = $out . "offset" . $offsets[$o] . "bp_" . $chrm . ".temp";
			$files{$chrm}{$offsets[$o]}{temp} = $cOut;
			my $pid = $pm->start and next;
			&process_chrm($chrm,$offsets[$o],$cOut,$count{$chr},$o);
			$pm->finish;
		}
	}

	if (!defined($alignDir)){
		$catCMD .= $cOut . " ";
		push(@delFiles,$cOut);
	}
}
$pm->wait_all_children;

if($alignDir){
	foreach my $chrm (keys %files){
		my $align_file = $alignDir . "ALIGN_" . $gb . "_" . $chrm . "_ADJUST.wig";
		foreach my $off (keys %{$files{$chrm}}){
			my $temp_input = $files{$chrm}{$off}{temp};
			$files{$chrm}{$off}{temp} .= ".align";
			$catCMD .= $files{$chrm}{$off}{temp} . " ";
			push(@delFiles,$files{$chrm}{$off}{temp});
			my $pid = $pm->start and next;
			&get_align($align_file,$temp_input,$files{$chrm}{$off}{temp});
			$pm->finish;
		}
	}
}
$pm->wait_all_children;

$catCMD .= " > $cnv_wins";
system($catCMD);
unlink($_) foreach @delFiles;

sub process_chrm{
	my ($chrm,$offset,$cOut,$count_ref,$o) = @_;
	my @hits = @{$count_ref->[0][$o]};
	my $numWins = @hits;
	open(OUT, ">$cOut");

	for my $id (0..$#hits){
		my $start = ($window_size*$id)+$offset+1;
		my $end = $start+$window_size-1;

		my $hits_raw = $hits[$id]+0;
		print OUT "$chrm\t$start\t$end\t$offset\t$hits_raw\n";
	}close OUT;
}

sub get_align {
	my ($alignFile,$tempWin,$winOut) = @_;
	open(ALIGN,$alignFile);
	my $currStart = 1;
	open(TEMP, $tempWin);
	open(OUT,">$winOut");
	
	while(<TEMP>){
		chomp;
		my $pLine = $_;
		my @line = split(/\t/,$_);
		if(!eof(ALIGN)){
			while($line[1] > $currStart){
				my $temp = <ALIGN>;
				$currStart++;
			}

			if($line[2] < $currStart){
				print OUT "$pLine\tNA\tNA\tstop before curr start\n";
			}elsif($line[1] == $currStart){
				my $alignCount = 0;
				while($currStart <= $line[2]){
					my $aScore = <ALIGN>;
					chomp($aScore);
					$currStart++;
					$alignCount += $aScore;
				}
				my $reads_alignbp = 0;
				$reads_alignbp = sprintf "%.5f", $line[4]/$alignCount if $alignCount > 0;
				print OUT "$pLine\t$alignCount\t$reads_alignbp\n";
			}
		}else{
			print OUT "$pLine\tNA\tNA\tend of align file\n";
		}
	}close TEMP;close OUT;
	close ALIGN;
	unlink($tempWin);
}

sub bst_readsTOcnv {
	my ($cnvStarts_aref,$cnvCoords_href,$readCoord_aref,$chrm,$tempOut) = @_;
	foreach my $read_coord (@{$readCoord_aref}){
		my $findOverlap = 0;
		my ($maxInt,$minInt) = ($#{$cnvStarts_aref},0);
		while( ($minInt <= $maxInt) && $findOverlap == 0){
			my $node = int(($maxInt + $minInt)/2);
			if($read_coord < ${$cnvStarts_aref}[$node]){
				$maxInt = $node - 1;
			}elsif($read_coord > $cnvCoords_href->{${$cnvStarts_aref}[$node]}{stop}){
				$minInt = $node + 1;
			}elsif( $read_coord >= ${$cnvStarts_aref}[$node] && $read_coord <= $cnvCoords_href->{${$cnvStarts_aref}[$node]}{stop}){
				$cnvCoords_href->{${$cnvStarts_aref}[$node]}{count}++;
				$findOverlap = 1;
			}
		}
	}
	
	open(OUT,">$tempOut");
	foreach my $start (@{$cnvStarts_aref}){
		my $stop = $cnvCoords_href->{$start}{stop};
		print OUT "$chrm\t$start\t$stop\t0\t" . $cnvCoords_href->{$start}{count} . "\n";
	}close OUT;
}
