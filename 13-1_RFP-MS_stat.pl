#!usr/bin/perl -w
use strict;

print << "EOF";
Screen circRNA support by both RFP and MS .
EOF
die "Usage: perl $0 <compendium.chr.sort.gene.exon.FLiso.RFP> <spectra.circ.filter> <OUT>\n" unless @ARGV==3;

open RFP, "<", "$ARGV[0]" or die ($!);
open MS, "<", "$ARGV[1]" or die ($!);
open OUT, ">", "$ARGV[2]" or die ($!);

my %MS;
my @RFP;

while(<RFP>){
    chomp;
    my @line = split(/\t/,$_);
    push @RFP,$line[3] if ($line[6]);
}
close RFP;

<MS>;
while(<MS>){
    chomp;
    my @line = split(/\t/,$_);
    my @proteins = split(/\//, $line[0]);
    foreach my $pro (@proteins){
        $pro =~ /(?<hg38>.+)\_(?<circID>hsa.+)\_(?<frame>.+)\:(?<start>\d+)\-(?<end>\d+)/;
        $MS{$+{circID}} = 1;
    }
}
close MS;

my @common = grep{$MS{$_}} @RFP;
@common = sort @common;
if(@common){
    print OUT join"\n",@common;
}else{
    print "no common circRNA support by both RFP and MS\n";
}
close OUT;

my @MS = keys %MS;
my $MS = @MS;
#print join"\t",@MS,"\n";
print "$MS circRNAs supported by MS\n";
