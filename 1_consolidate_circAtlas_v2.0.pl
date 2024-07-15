#usr/bin/perl -w
use strict;
print << "EOF";
This program is consolidate human circRNA information from circAtlas_v2.0.
EOF
die "Useage: perl $0 <humancircRNA.bed> <human_sequence.txt> <hg38_hg19.tsv> <compendium.tsv> \n" unless @ARGV==4;

open BED, "<","$ARGV[0]"or die ($!);
open SEQ, "<","$ARGV[1]"or die ($!);
open LINK, "<","$ARGV[2]"or die ($!);
open OUT, ">","$ARGV[3]"or die ($!);

my %circRNAs;

my $title_bed = <BED>;
chomp $title_bed;
while(<BED>){
    chomp;
    my @line = split(/\t/,$_);
    $circRNAs{$line[4]} = \@line;
}
close BED;
#my @circ= keys %circRNAs; my $circ = @circ; print $circ; #580718


my $blankcount;
my $title_link = <LINK>;
chomp $title_link;
while(<LINK>){
    chomp;
    my @line = split(/\t/,$_);
    if (exists $circRNAs{$line[1]}){
        for(my $i=0;$i<@line;$i++){
            push @{$circRNAs{$line[1]}},$line[$i];
        }
    }else{$blankcount += 1;}
}
close LINK;
print "not found link item: $blankcount\n";

while(<SEQ>){
    chomp;
    my @line = split(/\t/,$_);
    if (exists $circRNAs{$line[1]}){
        push @{$circRNAs{$line[1]}},$line[2];
    }
    else { print "no seq";push @{$circRNAs{$line[1]}},"-" ;}
}

close SEQ;

#print OUT $title_bed."\t"."mature_seq"."\t".$title_link."\n";
print OUT $title_bed."\t".$title_link."\t"."mature_seq"."\n";

my @circRNAs = keys %circRNAs;
foreach my $circRNA (@circRNAs){
    $circRNAs{$circRNA}[13] = "-" if (! $circRNAs{$circRNA}[13]); #the number of seqs in human_sequence.txt is less than circRNAs in humancircRNA.bed
#print "noseq" if ( ! $circRNAs{$circRNA}[13]);
    print OUT join"\t",@{$circRNAs{$circRNA}};
    print OUT "\n";
}
close OUT;
