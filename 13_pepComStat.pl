#!usr/bin/perl -w
use strict;
use List::Util qw(max);

print << "EOF";
This program is to get the maximum score of frame peptide for each circRNA.
EOF
die "Usage: perl $0 <class> <score> <OUT>\n" unless @ARGV==3;

open IN1, "<", "$ARGV[0]" or die ($!);
open IN2, "<", "$ARGV[1]" or die ($!);
open OUT, ">", "$ARGV[2]" or die ($!);

my %class;
while(<IN1>){
    chomp;
    my @line = split(/\t/,$_);
    if ($line[1] eq "positive"){
        $class{$line[0]} = 1;
    }
}
close IN1;

my %circ;
while(<IN2>){
    chomp;
    my @line = split(/\t/,$_);
    next unless (exists $class{$line[0]});
    my @info = split(/\,/,$line[0]);
    push @{$circ{$info[3]}}, $line[1];
}
close IN2;

my $count = 0;
foreach my $circ (keys %circ){
    my $max = max(@{$circ{$circ}});
    print OUT "$circ\t$max\n";
    $count += 1;
}
close OUT;
print "$count circRNAs(translated peptides) are positively similar to nature proteins\n";

