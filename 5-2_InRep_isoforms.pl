#!usr/bin/perl -w
use strict;
print << "EOF";

This program is to replace (and insert) circRNA isoforms with full-length exon coordinates in hg19(from xiaojuan), which have no exon been deleted by liftOver.

*Note: fullLength_exonHg38.bed & fullLength_exonHg38toHg19Unmap.bed are 1-based.(like gencode gtf)

EOF
die "Usage: perl $0 <compendium.chr.sort.gene.exon> <fullLength_exonHg38.bed> <fullLength_exonHg38toHg19Unmap.bed> <circRNA.exon.FLiso>\n" unless @ARGV==4;

open CIRC, "<", "$ARGV[0]" or die($!);
open HG38, "<", "$ARGV[1]" or die($!);
open UNMAP, "<", "$ARGV[2]" or die($!);
open OUT, ">", "$ARGV[3]" or die($!);

my %unmap;
while(<UNMAP>){
    chomp;
    next if($_ =~ /^\#/); #skip the discription part
    my @line = split(/\t/,$_);
    $unmap{$line[3]} = $line[4]; #key=isoformID, value=circID
}
close UNMAP;

my %isoform;
while(<HG38>){
    chomp;
    my @line = split(/\t/,$_);
    next if (exists $unmap{$line[3]}); #regard the isoform with deleted exons in hg38
    $isoform{$line[3]}[0] = $line[4]; #key=isoformID, value[0]=circID
    push @{$isoform{$line[3]}[1]},$line[1]; #@value[1]=exon_start
    push @{$isoform{$line[3]}[2]},$line[2]; #@value[1]=exon_end
}
close HG38;

my @isoforms = sort keys %isoform;
my %circRNA;
foreach my $isoform (@isoforms){
    my $starts = join"\,",@{$isoform{$isoform}[1]};
    my $ends = join"\,",@{$isoform{$isoform}[2]};
    my $circRNA = $isoform.":".$starts.":".$ends;
    push @{$circRNA{$isoform{$isoform}[0]}},$circRNA; #key=circID,value=isoformID_@starts_@ends;
}
#print join"\n",@{$circRNA{"hsa-RP11-206L10_0001"}};exit;

my $replace = 0;
my $repin_circ = 0;
my $repin_iso = 0;

while(<CIRC>){ #all circRNA here can cover isoforms' host circRNA
    chomp;
    my @line = split(/\t/,$_);
    my $circID = $line[3];
    my @newline;
    if (exists $circRNA{$circID}){
        if (@{$circRNA{$circID}} == 1){ #replace
            my @info = split(/\:/,${$circRNA{$circID}}[0]);
            my @start = split(/\,/,$info[1]);
            my @stop = split(/\,/,$info[2]);
            @start = sort { $a <=> $b } @start;
            @stop = sort { $a <=> $b } @stop;
            if ($start[0] >= $line[1] && $start[0] <= $line[2] && $start[-1] >= $line[1] && $start[-1] <= $line[2] && $stop[0] >= $line[1] && $stop[0] <= $line[2] && $stop[-1] >= $line[1] && $stop[-1] <= $line[2]){    
                my @temp = splice @line, 8, 5, ($info[1],$info[2],"-","-","-");
                print OUT join"\t",@line,"\n";
                $replace += 1;
            }else{
                print OUT join"\t",@line,"\n";    
            }
        }else { # @{$circRNA{$circID}} > 1, replace & insert
            my @starts; my @stops;
            for(my $i=0; $i<@{$circRNA{$circID}}; $i++){
                my @info = split(/\:/,${$circRNA{$circID}}[$i]);
                my @start = split(/\,/,$info[1]);
                my @stop = split(/\,/,$info[2]);
                push @starts,@start;
                push @stops,@stop;
            }
            @starts = sort { $a <=> $b } @starts;
            @stops = sort { $a <=> $b } @stops;

            if ($starts[0] >= $line[1] && $starts[0] <= $line[2] && $starts[-1]>= $line[1] && $starts[-1] <= $line[2] && $stops[0] >= $line[1] && $stops[0] <= $line[2] && $stops[-1] >= $line[1] && $stops[-1] <= $line[2]){
                $repin_circ += 1;
                for(my $i=0; $i<@{$circRNA{$circID}}; $i++){
                    my @info = split(/\:/,${$circRNA{$circID}}[$i]);
                    my @temp1 = splice @line, 3, 1, ($info[0]); #change $line[3], note keep $circID
                    my @temp2 = splice @line, 8, 5, ($info[1],$info[2],"-","-","-");
                    print OUT join"\t",@line,"\n";
                    $repin_iso += 1;
                }
            }else{
                print OUT join"\t",@line,"\n";
            }
        }
    }else{
        print OUT join"\t",@line,"\n";
    }
}
close CIRC;
close OUT;

my $iso = $repin_iso - $repin_circ;
print "$replace circRNA exon coordinates have been replaced\n";
print "$repin_circ circRNAs with $repin_iso isoforms have been inserted\n";
