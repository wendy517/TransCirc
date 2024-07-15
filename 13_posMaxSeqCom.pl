#!usr/bin/perl -w
use strict;
use List::Util qw(max);

print << "EOF";
This program is to get the maximum score of frame peptide for each circRNA.
EOF
die "Usage: perl $0 <score.class> <OUT>\n" unless @ARGV==2;

open IN, "<", "$ARGV[0]" or die ($!);
open OUT, ">", "$ARGV[1]" or die ($!);

my %posMaxSeqComp;
<IN>;
while(<IN>){
    chomp;
    my @line = split(/\t/,$_);
    if ($line[11] eq "positive"){
        push @{$posMaxSeqComp{$line[3]}}, $line[10]; #key=circID, @value=@score
    }
}
close IN;


my $count = 0;
print OUT "circID\tscore\n";

foreach my $circ (keys %posMaxSeqComp){
    my $max = max(@{$posMaxSeqComp{$circ}});
    print OUT "$circ\t$max\n";
    $count += 1;
}
close OUT;
print "$count circRNAs are positively similar to nature protein-coding sequences\n";

