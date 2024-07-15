#!usr/bin/perl -w
use strict;
print << "EOF";
This program is to score the association with ribosome/polysome of circRNA(source from independent study).
*Input bed format.
EOF

die "Usage: perl $0 <compendium.chr.sort.gene.exon.FLiso.seq.ORF.pep.row> <ZN2018NC_PMID_hg38> <XJF2017CR_PMID_hg19> <XJF2017CR_PMID_hg38> <ThomasPreiss2019SR_PMID_hg19> <ThomasPreiss2019SR_PMID_hg38> <hg38_hg19_v2.0.tsv> <compendium.chr.sort.gene.exon.FLiso.seq.ORF.pep.row.RFP>  \n" unless @ARGV==8;

open IN, "<", "$ARGV[0]" or die ($!);
open ZN38, "<", "$ARGV[1]" or die($!);
open XJF19, "<", "$ARGV[2]" or die($!);
open XJF38, "<", "$ARGV[3]" or die($!);
open TP19, "<", "$ARGV[4]" or die($!);
open TP38, "<", "$ARGV[5]" or die($!);
open LINK, "<", "$ARGV[6]" or die($!);
open OUT, ">", "$ARGV[7]" or die ($!);

my %ZN38;
my %XJF38;
my %XJF19;
my %TP19;
my %TP38;

sub circ_study{ #in: (STUDY); out:%STUDY
    my ($ref) = @_;
    my %ref;
    while(<$ref>){
        chomp;
        my @line = split(/\t/,$_);
        $ref{"$line[0]_$line[1]_$line[2]"} = 0; #key = chr_start_end
    }
    close $ref;
    return %ref;
}

%ZN38 = circ_study("ZN38");
%XJF38 = circ_study("XJF38");
%XJF19 = circ_study("XJF19");
%TP19 = circ_study("TP19");
%TP38 = circ_study("TP38");

my %link;
<LINK>; #skip title
while(<LINK>){
    chomp;
    my @line = split(/\t/,$_);
    $line[3] =~ /(?<chrom>chr.+)\:(?<start>\d+)\|(?<end>\d+)/;
    $link{$line[2]} = "$+{chrom}_$+{start}_$+{end}" if (exists $XJF19{"$+{chrom}_$+{start}_$+{end}"} || exists $TP19{"$+{chrom}_$+{start}_$+{end}"});
}
close LINK;

my $study = 3;
while(<IN>){
    chomp;
    my @line = split(/\t/,$_);
    next if($line[13] eq "-");
    my $count = 0;
    my @PMID = qw(-);   ################### update this in input requirment
    my @new_line;
    if ((exists $ZN38{"$line[0]_$line[1]_$line[2]"}) && ((exists $XJF38{"$line[0]_$line[1]_$line[2]"})||(exists $XJF19{$link{"$line[0]:$line[1]|$line[2]"}})) && ((exists $TP38{"$line[0]_$line[1]_$line[2]"})||(exists $TP19{$link{"$line[0]:$line[1]|$line[2]"}}))){
        $count = 3;
        @PMID = qw(30367041 26735365 28281539);
    }
    elsif( (exists $ZN38{"$line[0]_$line[1]_$line[2]"}) && ((exists $XJF38{"$line[0]_$line[1]_$line[2]"})||(exists $XJF19{$link{"$line[0]:$line[1]|$line[2]"}}))){
        $count = 2;
        @PMID = qw(30367041 28281539);
    }
    elsif( (exists $ZN38{"$line[0]_$line[1]_$line[2]"}) && ((exists $TP38{"$line[0]_$line[1]_$line[2]"})||(exists $TP19{$link{"$line[0]:$line[1]|$line[2]"}})) ){
        $count = 2;
        @PMID = qw(30367041 26735365);
    }
    elsif( ((exists $XJF38{"$line[0]_$line[1]_$line[2]"})||(exists $XJF19{$link{"$line[0]:$line[1]|$line[2]"}})) && ((exists $TP38{"$line[0]_$line[1]_$line[2]"})||(exists $TP19{$link{"$line[0]:$line[1]|$line[2]"}})) ){
        $count = 2;
        @PMID = qw(26735365 28281539);
    }
    elsif(exists $ZN38{"$line[0]_$line[1]_$line[2]"}){
        $count = 1;
        @PMID = qw(30367041);
    }
    elsif((exists $XJF38{"$line[0]_$line[1]_$line[2]"})||(exists $XJF19{$link{"$line[0]:$line[1]|$line[2]"}})){
        $count = 1;
        @PMID = qw(28281539);
    }
    elsif((exists $TP38{"$line[0]_$line[1]_$line[2]"})||(exists $TP19{$link{"$line[0]:$line[1]|$line[2]"}})){
    $count = 1;
    @PMID = qw(26735365);
    }
    @new_line = (@line[0..5], $count/$study, join"\,",@PMID);
    print OUT join"\t",@new_line,"\n";
}
close IN;
close OUT;
