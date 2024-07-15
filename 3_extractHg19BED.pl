#usr/bin/perl -w
use strict;
print << "EOF";
This program is to extract hg19 genomic position for each circRNA with customized input.
EOF
die "Usage: perl $0  <circRNA_hg19.input> <out.bed> \n" unless @ARGV==2;

open IN, "<","$ARGV[0]"or die ($!);
open OUT,">","$ARGV[1]"or die ($!);

while(<IN>){
    chomp;
    my @line = split(/\t/,$_);
    $line[13] =~ /(?<chrom>chr.+)\:(?<start>\d+)\|(?<end>\d+)/;
    print OUT "$+{chrom}\t$+{start}\t$+{end}\t$line[3]\t$line[7]\t$line[5]\n"; #chr start end circRNAID geneID strand
}
close IN;
close OUT;


