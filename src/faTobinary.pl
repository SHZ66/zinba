#!/usr/bin/perl

# faTobinary.pl - converts a fa file to 0 or 1s based on either GC or N
#
# faTobinary.pl listChrmFA.txt

if ($#ARGV != 1){
    print STDERR "\nEXITING: Incorrect number of parameters\nfaTobinary.pl listChrmFA.txt genomeBuild\n\n";
    exit;
}

my $listChrmFA = shift;
chomp($listChrmFA);

my $gb = shift;
chomp($gb);

my $clOut = $gb . "_c_chromLength.txt";
open(CLENGTH,">$clOut");

open(LIST,$listChrmFA);
while(<LIST>){
    chomp;
    my $gcout = $_;
#    my $nout = $_;
    $gcout =~ s/\..*/\.wig/g;
    $gcout = "GC_" . $gb . "_" . $gcout;
#    $nout =~ s/\..*/\.wig/g;
#    $nout = "N_" . $gb . "_" . $nout;
    open(FA,$_);
    open(GCOUT,">$gcout");
#    open(NOUT,">$nout");
    my $bpCount = 0;
    while(<FA>){
        chomp;
        if($_ =~ '>'){
            $_ =~ s/\>//g;
            print GCOUT "fixedStep chrom=$_ start=1 step=1\n";
#            print NOUT "fixedStep chrom=$_ start=1 step=1\n";
            print CLENGTH "$_\t1";
        }else{
            my $gcline = $_;
#            my $nline = $_;
            
            $gcline =~ s/[G|C]/1/g;
            $gcline =~ s/\D/0/g;
#            $nline =~ s/N/1/g;
#            $nline =~ s/\D/0/g;
            
            my @gcs = split("",$gcline);
#            my @ns = split("",$nline);
            $bpCount += ($#gcs + 1);
            
            print GCOUT join("\n",@gcs), "\n";
#            print NOUT join("\n",@ns), "\n";
        }
    }close FA;close GCOUT;#close NOUT;
    print CLENGTH "\t$bpCount\n";
}close LIST;close CLENGTH;