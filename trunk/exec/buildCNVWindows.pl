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
		--twoBit gb.2bit
		
		--window-size (default = 500 bp)
		--offset-size (default = 0)
		--perc-n-thresh (default 0.1)
		
#########################################

USAGE

my ($seq_hits,$cOut,%files);
my $window_size = 500;
my $offset = 0;
my $nThresh = 0.1;
my $twoBitFile = undef;
my $nProcesses = 1;

my $result = GetOptions(
	"seq=s" => \$seq_hits,
	"twoBit=s" => \$twoBitFile,
	'window-size=i' => \$window_size,
	"offset-size=i" => \$offset,
	"perc-n-thresh=f" => \$nThresh,
	"processes=i" => \$nProcesses,
);

die $usage unless($seq_hits);
my $pm = new Parallel::ForkManager($nProcesses);

my @offsets = 0;
my $numOffsets = 1;
$numOffsets = int($window_size/$offset) if $offset > 0;
for (my $o = 1; $o < $numOffsets; $o++){
	push(@offsets,($offset*$o));
}

open(SEQ,$seq_hits) or die;
my (%count,%chrom);
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

my $cnv_wins = $seq_hits;
my $out = $seq_hits;
$out =~ s/\..*//g;
$cnv_wins =~ s/\..*/\.cnvs/g;
my @delFiles;
my $catCMD = "cat ";

foreach my $chr(sort{$a<=>$b} keys %chrom){
	my $chrm = "chr" . $chr;
	for (my $o = 0; $o <= $#offsets; $o++){
		$cOut = $out . "offset" . $offsets[$o] . "bp_" . $chrm . ".temp";
		$files{$chr}{$offsets[$o]}{gcSeq} = $out . "offset" . $offsets[$o] . "bp_" . $chrm . ".gcseq";
		$files{$chr}{$offsets[$o]}{temp} = $cOut;
		if (!defined($twoBitFile)){
			$catCMD .= $cOut . " ";
			push(@delFiles,$cOut);
		}

                my $pid = $pm->start and next;
		&process_chrm($chrm,$offsets[$o],$cOut,$count{$chr},$o);
                $pm->finish;
	}
}
$pm->wait_all_children;

if(defined($twoBitFile)){
	foreach my $chr(sort{$a<=>$b} keys %chrom){
		my $chrm = "chr" . $chr;
		for (my $o = 0; $o <= $#offsets; $o++){
			my $gcSeq = $files{$chr}{$offsets[$o]}{gcSeq};
			my $start = $offsets[$o];
			$start-- if $start > 0;
			my $tempFile = $files{$chr}{$offsets[$o]}{temp};
			my $winOut = $tempFile . 2;
			$catCMD .= $winOut . " ";
			push(@delFiles,$winOut);

			my $cFileLen = $chrm . "_" . $offset . ".txt";
			`twoBitInfo /gbdb/hg18/hg18.2bit:$chrm $cFileLen`;
			open(CLEN,$cFileLen);
			my $line = <CLEN>;
			chomp($line);
			my ($cInfo,$cLength) = split(/\t/, $line);
			close CLEN; unlink($cFileLen);
			`twoBitToFa -noMask -seq=$chrm -start=$start -end=$cLength $twoBitFile $gcSeq 2> /dev/null`;

                        my $pid = $pm->start and next;
			&get_gcPerc($gcSeq,$tempFile,$winOut,$window_size);
                        $pm->finish;
		}
	}
	$pm->wait_all_children;
	

}

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
		print OUT "$chrm\t$start\t$end\t$hits_raw\n";
	}close OUT;
}

sub get_gcPerc {
	my ($gcFile,$tempWin,$winOut,$winSize) = @_;
	open(GC,$gcFile);
	my $gcHeader = <GC>;
	open(TEMP, $tempWin);
	open(OUT,">$winOut");
	my $readLen = 0;
	my $ncount = 0;
	while(<TEMP>){
		chomp;
                my $outLine = $_;
                my $gcFlag = 0;
                while($gcFlag == 0){
                        my $seq = <GC>;
                        chomp($seq);
                        $readLen += length($seq);
                        if($readLen == $winSize || eof(GC)){
                                $ncount += ($seq =~ tr/N//);
                                if($ncount/$winSize < $nThresh){
                                        print OUT "$outLine\n";
                                }
                                $readLen = 0;
                                $ncount = 0;
                                $gcFlag = 1;
                        }elsif($readLen < $winSize){
                                $ncount += ($seq =~ tr/N//);
                        }elsif($readLen > $winSize){
                                my $tempSeq = substr($seq,($readLen-$winSize));
                                my $seq = substr($seq,0,($readLen-$winSize-1));
                                $ncount += ($seq =~ tr/N//);
                                if($ncount/$winSize < $nThresh){
                                        print OUT "$outLine\n";
                                }
                                $readLen = length($tempSeq);
                                $ncount = ($tempSeq =~ tr/N//);
                                $gcFlag = 1;
                        }
	}close TEMP;close OUT;
	close GC;
	unlink($tempWin);unlink($gcFile);
}
