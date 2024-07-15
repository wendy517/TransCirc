#usr/bin/perl -w
use strict;
print << "EOF";
This program is to get gene.bed6(with gene_name, gene_id and strand) from gencode.gtf.
EOF
die "Usage: perl $0 <gencode.gtf> <gene.bed6> \n" unless @ARGV==2;

open GTF,"<","$ARGV[0]"or die ($!);
open BED,">","$ARGV[1]"or die ($!);

my %genes;
while(<GTF>){
    chomp;
    next if($_ =~/^\#/); #skip description part of gtf
    my @line = split(/\t/,$_);
    if ($line[2] eq "gene"){
        $line[8] =~/gene_id \"(?<gene_id>ENSG.+)\"\; gene_type \"(.+)\"\; gene_name \"(?<gene_name>.+)\"\; level (\d+)\;.*/; #for mouse:ENSMUSG; for human:ENSG
        $genes{$+{gene_id}}[0] = $line[0]; #chr
        $genes{$+{gene_id}}[1] = $line[3]; #start
        $genes{$+{gene_id}}[2] = $line[4]; #end
        $genes{$+{gene_id}}[3] = $+{gene_name};
        $genes{$+{gene_id}}[4] = $line[6]; #strand
    }
}
close GTF;

my @gene_ids = keys %genes;
foreach my $gene_id (@gene_ids){
    print BED "$genes{$gene_id}[0]\t$genes{$gene_id}[1]\t$genes{$gene_id}[2]\t$genes{$gene_id}[3]\t$gene_id\t$genes{$gene_id}[4]\n";    
}
close BED;
