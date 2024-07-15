#!usr/bin/perl -w
use strict;
print << "EOF";
This program is to score circRNAs according to IRES hexmer(s) searched and mergerd in circRNAs sequence. 
EOF
die "Usage: perl $0 <circRNA.seq> <IRES.csv> <circRNA.seq.IRES>\n" unless @ARGV==3;

open CIRC, "<", "$ARGV[0]" or die ($!);
open IRES, "<", "$ARGV[1]" or die ($!);
open OUT, ">", "$ARGV[2]" or die ($!);

#######store IRES hexmer score
my %Zhex;
while(<IRES>){
    chomp;
    my @line = split(/\,/,$_);
    $Zhex{$line[0]} = $line[6]; #key = hexmer, value = Zscore
}
close IRES;

my @IRESs = sort {$Zhex{$b} <=> $Zhex{$a}} keys %Zhex; #sort by Zscore, decreasing = T
#print join"\t",@IRES,"\n";exit;

#######search and merge
my $count = 0;
while(<CIRC>){
    chomp;
    my @line = split(/\t/,$_);
    next if ($line[13] eq "-");
    my @nt = split(//,$line[13]);
    
    print OUT join"\t",(@line[0..5], 0, "-", "-"),"\n" if (@nt < 6); #length($seq < 6)
    next if (@nt < 6);

    my $seq = $line[13].$nt[0].$nt[1].$nt[2].$nt[3].$nt[4] ; #across BSJ
    my %start;
    my @seq_hex;
    for(my $i=0; $i<length($seq)-5; $i++){
        my $hex = substr($seq, $i, 6);
        push @seq_hex, $hex;
    }
    foreach my $IRES (@IRESs){
        my @index = grep $seq_hex[$_] eq $IRES, 0..$#seq_hex;
        foreach my $index (@index){
            $start{$index} = $Zhex{$IRES}; #key=index_hex, value = score_hex
        }
    }

    my @starts = sort {$a <=> $b} keys %start;
    print OUT join"\t",(@line[0..5], 0, "-", "-"),"\n" if (!@starts); #empty IRES index
    next unless (@starts);

    ###@starts; @ends; @total_scores;
    my @ends; my @scores;
    for(my $i=0; $i<@starts; $i++){
        $ends[$i] = $starts[$i] + 5;
#        $ends[$i] -= length($line[13]) if ($ends[$i]+1 > length($line[13]));
        $scores[$i] = $start{$starts[$i]};
    }
#print join"\t", @starts ,"\n";exit;
    ###merge
    foreach my $i (0..$#starts-1){
        my $j = $i + 1;
        if($ends[$i] >= $starts[$j]){
            $starts[$j] = $starts[$i];
            $scores[$j] += $scores[$i];
            $starts[$i] = $ends[$i] = $scores[$i] = 0;
        }
    }
    ###calculate total IRES score
    my $IRESscore = 0;
    my @start; my @end;
    foreach my $i (0..$#scores){
        next if ($scores[$i] == 0);
        my $score = $scores[$i]/($ends[$i] - $starts[$i] + 1);
        $IRESscore += $score;
        push @start,$starts[$i];
        $ends[$i] -= length($line[13]) if ($ends[$i]+1 > length($line[13]));
        push @end,$ends[$i];
    }
    my $starts = join"\,", @start;
    my $ends = join"\,", @end;
    print OUT join"\t",(@line[0..5], $IRESscore, $starts, $ends),"\n";
    $count += 1;
}
close CIRC;
close OUT;
print "$count circRNA(s) have IRES score record\n" ;

################log10
sub log10 {
    my $n = shift;
    return log($n)/log(10);
}
