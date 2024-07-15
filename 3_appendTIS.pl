#usr/bin/perl -w
use strict;
print << "EOF";
This program is to append TIS information(like exon coordinates in bed) for each circRNAs with intersected file.
EOF
die "Usage: perl $0  <input> <intersect> <out> \n" unless @ARGV==3;

open IN, "<","$ARGV[0]"or die ($!);
open INT, "<","$ARGV[1]"or die ($!);
open OUT, ">","$ARGV[2]"or die ($!);

my %TISs; #key= circID_geneID
while(<INT>){
    chomp;
    my @line = split(/\t/,$_);
    push @{$TISs{$line[9]."_".$line[10]}[0]}, $line[1]; #Coordinate start codon
    push @{$TISs{$line[9]."_".$line[10]}[1]}, $line[2]; #Coordinate stop codon
    push @{$TISs{$line[9]."_".$line[10]}[2]}, $line[3]; #TIS type
    push @{$TISs{$line[9]."_".$line[10]}[3]}, $line[4]; #Predicted uORF
}
close INT;

while(<IN>){
    chomp;
    my @line = split(/\t/,$_);
    if (exists $TISs{$line[3]."_".$line[7]}){
        my $Coordinate_start_codon = join"\,",@{$TISs{$line[3]."_".$line[7]}[0]};
        my $Coordinate_stop_codon = join"\,",@{$TISs{$line[3]."_".$line[7]}[1]};
        my $TIS_type = join"\,",@{$TISs{$line[3]."_".$line[7]}[2]};
        my $Predicted_uORF = join"\,",@{$TISs{$line[3]."_".$line[7]}[3]};
        my @temp = splice @line, 10, 0, ($Coordinate_start_codon,$Coordinate_stop_codon,$TIS_type,$Predicted_uORF); #insert before parameter2
    }else {my @temp = splice @line, 10, 0, qw(- - - -) ;}
    print OUT join"\t",@line,"\n";
}
close IN;
close OUT;
