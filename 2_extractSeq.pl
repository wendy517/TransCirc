#!usr/bin/perl -w
use strict;
print << "EOF";
This program is to collect exon corrdinates of circRNAs and extract the correspoding sequence from genome (bedtools getfasta && cat)
*CircRNAs with partial seq or no exon are excluded.
EOF
die "Usage: perl $0 <cirRNA.exon> <genome> <out.seq> \n" unless @ARGV==3;

open EXON,"<","$ARGV[0]"or die ($!);
open HG38,"<","$ARGV[1]"or die ($!);
open OUT,">","$ARGV[2]"or die ($!);

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
#####################append sequence
my $nothing = "-";
while(<EXON>){
    chomp;
    my @line = split(/\t/,$_);
    my @starts = split(/\,/,$line[8]); #exon_starts
    my @ends = split(/\,/,$line[9]); #exon_ends
#print join "\t",@starts,"\n";
#print join "\t",@ends,"\n";
    my $seq;
    if ($nothing!~@starts){ #exon != "-"
        if ($line[5] eq "+"){
            for (my $i=0; $i<@starts; $i++){
                my $start = $starts[$i];
                my $end = $ends[$i];
                my $length = $end - $start + 1;
#print $length,"\n";
                my $exon_seq = substr($chr{$line[0]}, $start-1, $length);
#print $exon_seq ,"\n";
                $seq .= $exon_seq;    
            }
#print $seq,"\n";
        }
        if ($line[5] eq "-"){ #re-order the @exon for neg strand, then reverse
            for (my $i=$#starts; $i>=0; $i--){
                my $start = $starts[$i];
                my $end = $ends[$i];
                my $length = $end - $start + 1;
                my $exon_seq = substr($chr{$line[0]}, $start-1, $length);
                $seq .= $exon_seq;
            }
                $seq  =~tr/ACGTacgt/TGCAtgca/;
                $seq = reverse($seq);
        }
        
        if ($seq){push @line, $seq;}
        else{push @line, "-";}
#print $seq,"\n";exit;
    }
    else {push @line, "-";}
    print OUT join"\t",@line,"\n";
}
close OUT;
