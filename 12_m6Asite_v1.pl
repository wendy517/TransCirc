#!usr/bin/perl -w
use strict;
print << "EOF";
This program is to append m6A site(s) whoese peaks detected by 3 tools (exomePeak, MACS2 and MeTPeak) in circRNA sequence(exon).
EOF

die "Usage: perl $0 <circRNA.seq> <m6A.merge> <circRNA.m6A> <hg38.fa>\n" unless @ARGV==4;

open CIRC, "<", "$ARGV[0]" or die ($!);
open M6A, "<", "$ARGV[1]" or die ($!);
open OUT, ">", "$ARGV[2]" or die ($!);
open HG38, "<", "$ARGV[3]" or die ($!);

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

#################store m6A data from REPIC
my %mergePeak;
while(<M6A>){
    chomp;
    my @line = split(/\t/,$_);
    push @{$mergePeak{$line[0]}[0]}, $line[1]; #key = chr, value[0] = start
    push @{$mergePeak{$line[0]}[1]}, $line[2]; #key = chr, value[1] = end
}
close M6A;


######################A in peak within exon
my $count = 0;
my $num_m6A;
while(<CIRC>){
    chomp;
    my @line = split(/\t/,$_);
    next if ($line[13] eq "-");

    my @exonStart = split(/\,/,$line[8]);
    my @exonEnd = split(/\,/,$line[9]);
    my @exonPeakStart;
    my @exonPeakEnd;
    for(my $i=0; $i<@exonStart; $i++){
        for(my $j=0; $j<@{$mergePeak{$line[0]}[0]}; $j++){
            next if ($exonStart[$i] > ${$mergePeak{$line[0]}[0]}[$j] || $exonEnd[$i] < ${$mergePeak{$line[0]}[1]}[$j]); #start1>end2, end1<start2, no overlap
            if ($exonStart[$i] < ${$mergePeak{$line[0]}[0]}[$j] && $exonEnd[$i] > ${$mergePeak{$line[0]}[1]}[$j]){ #1 cover 2
                my $start = ${$mergePeak{$line[0]}[0]}[$j];
                my $end = ${$mergePeak{$line[0]}[1]}[$j];
                push @exonPeakStart, $start;
                push @exonPeakEnd, $end;
            }
            elsif($exonStart[$i] > ${$mergePeak{$line[0]}[0]}[$j] && $exonEnd[$i] < ${$mergePeak{$line[0]}[1]}[$j]){  #2 cover 1
                my $start = $exonStart[$i];
                my $end = $exonEnd[$i];
                push @exonPeakStart, $start;
                push @exonPeakEnd, $end;
            }
            elsif($exonStart[$i] > ${$mergePeak{$line[0]}[0]}[$j] && $exonEnd[$i] > ${$mergePeak{$line[0]}[1]}[$j]){ #2<1
                my $start = $exonStart[$i];
                my $end = ${$mergePeak{$line[0]}[1]}[$j]; 
                push @exonPeakStart, $start;
                push @exonPeakEnd, $end;
            }
            elsif($exonStart[$i] < ${$mergePeak{$line[0]}[0]}[$j] && $exonEnd[$i] < ${$mergePeak{$line[0]}[1]}[$j]){ #1<2
                my $start = ${$mergePeak{$line[0]}[0]}[$j];
                my $end = $exonEnd[$i];
                push @exonPeakStart, $start;
                push @exonPeakEnd, $end;
            }
        }
    }

    print OUT join"\t",(@line[0..5],"-"),"\n" if (!(@exonPeakStart)); #no peak in exon
    next if (!(@exonPeakStart));

    my @Asites;
    for (my $i=0; $i<@exonPeakStart; $i++) {
        my $length = $exonPeakEnd[$i] - $exonPeakStart[$i] + 1;
        my $str = substr($chr{$line[0]}, $exonPeakStart[$i]-1, $length);
#        if ($line[5] eq "-"){
#            $str =~ tr/ACGTacgt/TGCAtgca/;
#            $str = reverse($str);
#        }
        my @tri_nt ;
        my @index;
        for(my $i=0; $i<length($str)-3; $i++){
            my $tri = substr($str, $i, 3);
            push @tri_nt, $tri;
        }
#        @index = grep $tri_nt[$_] =~/(A|G)AC/, 0..$#tri_nt; #grep all RAC
        @index = grep $tri_nt[$_] =~/(A|G)AC/, 0..$#tri_nt if ($line[5] eq "+"); #grep all RAC
        @index = grep $tri_nt[$_] =~/GT(T|C)/, 0..$#tri_nt if ($line[5] eq "-"); #grep all RAC in - strand
        next if (!@index);
        foreach my $index (@index){
            push @Asites, $index  + $exonPeakStart[$i] + 1; #tri+1, fif+2
            $num_m6A += 1;
        }
    }

    print OUT join"\t",(@line[0..5],"-"),"\n" if (!@Asites); #no A in peak
    next if (!@Asites);

    my $site = join"\,",@Asites;
    print OUT join"\t",(@line[0..5], $site),"\n";
    $count += 1;
}
close OUT;
print "$count circRNA(s) have $num_m6A m6A site(s) records\n";
