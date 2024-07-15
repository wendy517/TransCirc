#!usr/bin/perl -w
use strict;
print << "EOF";
This program is to interpolate genomic information extracting from query bed and gencode gtf into compendium.tsv.
EOF
die "Usage: perl $0 <compendium_queryBed.tsv> <gecode.gtf> <out.tsv> <Species:hsa>\n" unless @ARGV==4;

open GTF,"<","$ARGV[1]"or die ($!);
open BED,"<","$ARGV[0]"or die ($!);
open OUT,">","$ATGV[2]"or die ($!);

my %host_genes;
while(<BED>){
    chomp;
    my @line = split(/\t/,$_);
    $line[4] =~ /$AGRV[3]\-(?<gene_Name>)\_(\d+)/;
    $host_genes{$+{gene_Name}} = 1; 
}
close BED;

my %genes;
my @exons;
my @exonStarts;
my @exonEnds;

while(<GTF>){
    chomp;
    next if($_ =~/^\#/); #skip description part of gtf
    my @line = split(/\t/,$_);
    if ($line[2] eq "gene"){
        $line[8] =~/gene_id \"(?<gene_id>ENSG.+)\"\; gene_type \"(.+)\"\; gene_name \"(?<gene_name>.+)\"\; level (\d+)\;.*/; #for mouse:ENSMUSG; for human:ENSG
        if (exists $host_genes{$+{gene_name}}){ #filter genes
            $genes{$+{gene_name}."_".$+{gene_id}}[0] = $line[0]; #chr
            $genes{$+{gene_name}."_".$+{gene_id}}[1] = $line[3]; #start
            $genes{$+{gene_name}."_".$+{gene_id}}[2] = $line[4]; #end
            $genes{$+{gene_name}."_".$+{gene_id}}[3] = $+{gene_name};
            $genes{$+{gene_name}."_".$+{gene_id}}[4] = $+{gene_id};
        }
    }
    if ($line[2] eq "exon"){    
        $line[8] =~ /gene_id \"(?<geneID>ENSG.+)\"\; transcript_id \"(.+)\"\; gene_type \"(.+)\"\; gene_status \"(.+)\"\; gene_name \"(?<geneName>.+)\"\; .*/;
        if (exists $genes{$+{geneName}."_".$+{geneID}}){
            my $exon_coordinate = $line[0]."_".$line[3]."_".$line[4].$+{geneName}."_".$+{geneID}; #exon:chr_start_end_geneName_geneID
            push @exons,$exon_coordinate;
        }
    }
}
close GTF;
#my @genes_ids = keys %genes; #geneName_id

my %exons_gene;
foreach my $exon (@exons){
    my @info = split(/\_/,$exon);
    push @{$exons_gene{$info[3]."_".$info[4]}},$info[1]."_".$info[2]; #group exon pairs by geneName_id   
}

##################### join sites pair #####################
sub joinPair{
    my ($ref_s,$ref_e) = @_;
    my @start = @$ref_s;
    my @end = @$ref_e;
    my @pairs;
    if ($#start == $#end){
        for (my $i=0;$i<@end;$i++){
            my $pair = "$start[$i]_$end[$i]";
            push @pairs,$pair;
        }
    }else{print "not pairs!\n";}
    return @pairs;
}
##################### split sites pair #####################
sub splitPair{
    my ($ref) = @_;
    my @pairs = @$ref;
    my @starts;
    my @ends;
    for (my $i=0;$i<@pairs;$i++){
        my @pair = split("_",$pairs[$i]);
        push @starts,$pair[0];
        push @ends,$pair[1];
    }
    return (\@starts,\@ends);
}
##################### merge overlapped region #####################
sub mergeRegion{ #only accept array sorted by start site
    my ($ref_s,$ref_e) = @_;
    my @start = @$ref_s;
    my @end = @$ref_e;
    if ($#end == $#start && $#end > 0 ){
        foreach my $i (0..$#end - 1){
            my $j = $i + 1;
            if ($end[$i] >= $end[$j]){
                $end[$j] = $end[$i];
                $start[$j] = $start[$i];
                $start[$i] = $end[$i] = 0;
            }
            elsif($end[$i] >= $start[$j]){
                $start[$j] = $start[$i];
                $start[$i] = $end[$i] = 0;
            }
        }
    }
    @starts = grep(/^[1-9]/, @start); #remove 0
    @ends = grep(/^[1-9]/, @end); #remove 0
    return (\@starts,\@ends);
}

########################### merge exons per gene in gtf #########################
my %uniqexons_gene;
my @gene_ids = keys %exons_gene;
foreach my $gene_id (@gene_ids){
    my %count;
    my @{$uniqexons_gene{$gene_id}} = grep { ++$count{$_} == 1 } @{$exons_gene{$gene_id}}; #remove duplicate exons per gene
    @{$uniqexons_gene{$gene_id}} = map {$_ ->[0]} sort {$a->[1] <=> $b->[1]} map {[$_,$_=~/(\d+)_\d+/]} @{$uniqexons_gene{$gene_id}}; #sort by start site
    my ($ref_eSs,$ref_eEs) = splitPair(\@{$uniqexons_gene{$gene_id}});
    my @exonStarts = @$ref_eSs;
    my @exonEnds = @$ref_eEs;
    ($ref_eSs,$ref_eEs) = mergeRegion(\@exonStarts,\@exonEnds); #merge exons
    @exonStarts = @$ref_eSs;
    @exonEnds = @$ref_eEs;
    @{$uniqexons_gene{$gene_id}} = joinPair(\@exonStarts,\@exonEnds);
}


########################query##############################
my @genes_ids = keys %genes;
while(<BED>){
    chomp;
    my @line = split(/\t/,$_);
    $line[4] =~ /$AGRV[3]\-(?<gene_Name>)\_(\d+)/; 
    my $key = 
    if (exists )
        
}
