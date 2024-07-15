#!usr/bin/perl -w
use strict;

print << "EOF";
This program is to combine human-originated(OX=9606) peptide sequence(.fa) with potential peptide translated by circRNA.
60c in one row.
EOF

die "Usage: perl $0 <uniprot_sprot.fa> <circ.pep.column> <cmb.fa> \n" unless @ARGV==3;

open UNIP, "<", "$ARGV[0]" or die ($!);
open CIRC, "<", "$ARGV[1]" or die ($!);
open OUT, ">", "$ARGV[2]" or die ($!);

$/ = ">";
while(<UNIP>){
    chomp;
    next unless ($_ =~ /OX=9606/);
    print OUT "\>$_";
}
close UNIP;

$/ = "\n";
while(<CIRC>){
    chomp;
    my @line = split(/\t/,$_);
    next if($line[11] eq "-");
    my $title = "$line[0]\:$line[1]\-$line[2]\($line[5]\)\_$line[3]\_$line[6]\:$line[7]\-$line[8]"; #chr:start-end(strand)_circID@frame:orfStart-orfEnd
    print OUT "\>$title","\n";
    $line[11] =~ s/\*//g; #remove the stop codon
    for (my $i = 0; $i < length($line[11]); $i+=60){
        print OUT substr($line[11], $i, 60),"\n";
    }
}
close CIRC;
close OUT;
