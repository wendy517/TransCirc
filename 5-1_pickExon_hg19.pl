#!usr/bin/perl -w
use strict;
print << "EOF";
This program is to pick out the exon coordinates in hg19(bed6) for overlapped circRNA with full-length support that can be processed through liftOver.
*CircRNA isoforms are sorted by count. (decreasing = T)
Example: perl 5-1_pickExon_hg19.pl compendium.chr.sort.gene hg38_hg19_v2.0.tsv merge.AS.ciri.anno.bed fullLength_exonHg19.bed (1-based, which consistant with coordinates in gtf)
EOF
die "Usage: perl $0 <circRNA> <hg38_hg19_v2.0.tsv> <merge.AS.ciri.anno.bed> <exon_hg19.bed> <lost_circRNA_hg19>\n" unless @ARGV==5;

open CIRC, "<", "$ARGV[0]" or die($!);
open LINK, "<", "$ARGV[1]" or die($!);
open FLhg19, "<", "$ARGV[2]" or die($!);
open OUT, ">", "$ARGV[3]" or die($!);
open LOST, ">", "$ARGV[4]" or die($!);

my %circRNA;
while(<CIRC>){
    chomp;
    my @line = split(/\t/,$_);
    $circRNA{$line[3]} = $line[5]; #key = circAltas_ID, value = strand 408882
}
close CIRC;

my $seq = 0;
<LINK>;
my %BSJ_hg19;
while(<LINK>){
    chomp;
    my @line = split(/\t/,$_);
    next unless (exists $circRNA{$line[1]}); 
    $seq += 1;
    $BSJ_hg19{$line[3]}[0] = $line[1]; #key=BSJ_hg19, value[0] = circAltas_ID
    @{$BSJ_hg19{$line[3]}[1]} = ();
}
close LINK;

my @lost;
my $lost_isoform = 0;
while(<FLhg19>){
    chomp;
    my @line = split(/\t/,$_);
    next if ($line[7] eq "ICF" || $line[9]==0); # make sure intron-free
    my $bsj_hg19 = $line[0].":".($line[1]+1)."|".$line[2];
    if (exists $BSJ_hg19{$bsj_hg19}){
        my $isoform = $line[4]."_".$bsj_hg19."_".$line[9]."_".$line[10]."_".$line[11]; 
        push @{$BSJ_hg19{$bsj_hg19}[1]},$isoform; #value[1]= count_bsjHg19_exonCount_lengths_relativeStarts
    }
    else{
        $lost_isoform +=1;
        next if grep {$bsj_hg19 eq $_ } @lost;
        push @lost,$bsj_hg19;
    };
}
close FLhg19;
print LOST join"\n",@lost;
#foreach my $i (@{$BSJ_hg19{"chr17:556531|602620"}[1]}){
#    print $i,"\n";
#}exit;

my @BSJ_hg19 = keys %BSJ_hg19;
my $circRNAs = 0;
##########################################################sorting
foreach my $bsj_hg19 (@BSJ_hg19){
#print $bsj_hg19,"\n";exit;
    if (@{$BSJ_hg19{$bsj_hg19}[1]}){
#print $bsj_hg19,"\n";exit;
        $circRNAs += 1;
        @{$BSJ_hg19{$bsj_hg19}[1]} = map {$_ ->[0]} sort {$b->[1] <=> $a->[1]} map {[$_,$_=~/(\d+)\_.+/]} @{$BSJ_hg19{$bsj_hg19}[1]}; #sort by CIRI count callout, then give them serial number in order
        foreach my $circ_isoform (@{$BSJ_hg19{$bsj_hg19}[1]}){ #yes this is in order
            my @circ = split(/\_/,$circ_isoform);
            my $coord = $circ[2]."_".$circ[3]."_".$circ[4];
            push @{$BSJ_hg19{$bsj_hg19}[2]},$coord;
        }
    }
}
#########################################################remove duplicate isoform and print OUT
my $isoforms = 0;
foreach my $bsj_hg19 (@BSJ_hg19){
    my %count;
    @{$BSJ_hg19{$bsj_hg19}[2]} = grep { ++$count{$_} < 2 } @{$BSJ_hg19{$bsj_hg19}[2]}; #remove duplicated isoforms
    $isoforms += ($#{$BSJ_hg19{$bsj_hg19}[2]}+1);
    $bsj_hg19 =~ /(?<chr>chr.+)\:(?<start_plus1>\d+)\|(?<end>\d+)/;
    for (my $i=0; $i<@{$BSJ_hg19{$bsj_hg19}[2]}; $i++){ #yes this is in orde
        my @coord = split(/\_/,${$BSJ_hg19{$bsj_hg19}[2]}[$i]);
        my @starts = split(/\,/,$coord[2]);
        my @length = split(/\,/,$coord[1]);
        for(my $j=0; $j<$coord[0]; $j++){
            my $start = $starts[$j] + $+{start_plus1};
            my $end = $start + $length[$j];
            my $name = $BSJ_hg19{$bsj_hg19}[0]."-".$i;
            print OUT "$+{chr}\t$start\t$end\t$name\t$BSJ_hg19{$bsj_hg19}[0]\t$circRNA{$BSJ_hg19{$bsj_hg19}[0]}\n"; #chr,start,end,isoName,circName,strand
        }

    }
}

print "$seq circRNAs with full-length sequence in circAtlas\n"; #408882
print "$circRNAs circRNA with $isoforms isoforms with full-length need to replace sequences\n";
print $#lost+1," lost circRNAs(may locate in intergenic) with $lost_isoform lost isoforms with full-length no need to replace\n";
#sort -k1,1 -k2,2n
