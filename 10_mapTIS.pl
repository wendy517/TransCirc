#!usr/bin/perl -w
use strict;
print << "EOF";
This program is to map TIS data from TISdb to circRNAs.
EOF
die "Usage: perl $0 <circRNA> <hg38genome> <TIS_hg38lift> <circRNA.TIS> \n" unless @ARGV==4;

open CIRC, "<", "$ARGV[0]" or die ($!);
open HG38, "<", "$ARGV[1]" or die ($!);
open TIS, "<", "$ARGV[2]" or die ($!);
open OUT, ">", "$ARGV[3]" or die ($!);

#################store genome sequence
my %chr;
my $chrom;
while(<HG38>){
    chomp;
    if ($_ =~ /^>(..*)/){
        $chrom = $1;
    }
    else{
        $chr{$chrom} .= $_;
    }
}
close HG38;

################select TIS
my %TIS;
while(<TIS>){
    chomp;
    my @line = split(/\t/,$_);
    next if ($line[0] =~ /\_/); #scaffold in chrom
    my $codon;
    if ($line[5] eq "+"){
        $codon = substr($chr{$line[0]}, $line[1]-1, 3);
    }
    elsif($line[5] eq "-"){
        $codon = substr($chr{$line[0]}, $line[1]-1, 3);
        $codon =~ tr/ACGTacgt/TGCAtgca/;
        $codon = reverse($codon);
    }
    next if ($codon ne $line[3]); #ignore liftover mistakes
    push @{$TIS{$line[0]}},$line[1]; #@{$TIS{chrX}} = (TIS1, TIS2, ...)
}
close TIS;

################append TIS within circRNAs'(exon) range
my $count = 0;
my $num_TIS=0;
while(<CIRC>){
    chomp;
    my @line = split(/\t/,$_);
    next if ($line[13] eq "-");
    my @TISs = ();
    my @new_line;
    my @exonStart = split(/\,/,$line[8]);
    my @exonEnd = split(/\,/,$line[9]);
    foreach my $start (@{$TIS{$line[0]}}){
        for(my $i=0; $i<@exonEnd; $i++){
            if ($start >= $exonStart[$i] && $start+2 <= $exonEnd[$i]){
                push @TISs, $start;
                $num_TIS += 1;
            }
        }
    }
    if (@TISs){
        $count += 1;
        my $TISs = join"\,",@TISs;
        @new_line = (@line[0..5],$TISs);
    }else{
        @new_line = (@line[0..5],"-");
    }
    print OUT join"\t",@new_line,"\n";
}
close CIRC;
close OUT;
print "$count circRNAs have $num_TIS anotated alternative TIS in TISdb\n";
