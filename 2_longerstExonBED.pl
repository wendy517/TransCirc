#!usr/bin/perl -w
use strict;
print << "EOF";
This program is to extract and append the exon coordinates within circRNAs from the longest transcript in gencode.gtf for the input file with gene_id.
EOF
die "Usage: perl $0 <gencode.gtf> <input.gene> <out.gene.exon> \n" unless @ARGV==3;

open GTF,"<","$ARGV[0]"or die ($!);
open GENE,"<","$ARGV[1]"or die ($!);
open OUT,">","$ARGV[2]"or die ($!);

my %genes;
my %circRNA;
while(<GENE>){
    chomp;
    my @line = split(/\t/,$_);
    my @gene_id = split(/\,/,$line[7]);
    for(my $i = 0; $i < @gene_id; $i++){
        $genes{$gene_id[$i]} = (); #key=gene_ID
        push @{$circRNA{$line[3]}[0]},$gene_id[$i]; #key = circRNAID, value[0] = gene_id
    }
}
close GENE;
my @circRNAs = keys %circRNA;
#print $#circRNAs+1; exit;
#print join"\,",@{$circRNA{"hsa-A1BG-AS1_0001"}[0]};exit;
#######################get all tids and its length
while(<GTF>){
    chomp;
    next if($_ =~/^\#/); #skip description part of gtf
    my @line = split(/\t/,$_);
    if($line[2] eq "transcript"){
        $line[8] =~ /gene_id \"(?<gene_id>ENSG.+)\"\; transcript_id \"(?<transcript_id>ENST.+)\"\; gene_type \"(.+)\"\; gene_status \"(.+)\"\; gene_name \"(?<geneName>.+)\"\; transcript_type \"(.+)\"\; transcript_name \"(.+)\"\; .*/;
        next unless (exists $genes{$+{gene_id}});
        my $length = $line[4]-$line[3]; # open 1-based for gencode gtf: (1,5)
        my $length_transcriptID = $length."_".$+{transcript_id};
#        print $length_transcriptID ; exit;
        push @{$genes{$+{gene_id}}},$length_transcriptID; #key=gene_ID, value = all length_transcriptID
    }     
}
close GTF;
#my @keys = keys %genes;
#print join"\,",@keys;exit;
#foreach my $key (@keys){
#    print join"\,",@{$genes{$key}},"\n";
#}
#exit;
#######################pick the longest tid for each gene
my @gene_ids = keys %genes;
my %longest_tid;
foreach my $gene_ID (@gene_ids){
    @{$genes{$gene_ID}} = map {$_ ->[0]} sort {$b->[1] <=> $a->[1]} map {[$_,$_=~/(\d+)\_.+/]} @{$genes{$gene_ID}}; #sort by length, decreasing = T
    my @temp = split(/\_/,$genes{$gene_ID}[0]);
#    print join"\,",@temp,"\n";exit;
    $longest_tid{$gene_ID} = $temp[0]."_".$temp[1] ; # key = geneID, value = length_tid
#    print $longest_tid{$gene_ID},"\n";exit;
}    
#my @keys = keys %longest_tid; 
#print join"\,",@keys;exit;
#foreach my $key (@keys){
#    print "$key\t$longest_tid{$key}\n";
#}
#exit;
#########################pick transcript and according gene for circRNAs with multi-gene annotations
my @circRNA = keys %circRNA;
my @filtered_gene_ids;
foreach my $circRNA (@circRNA){
#$circRNA{$circRNA}[1] = ${$circRNA{$circRNA}[0]}[0] if(@{$circRNA{$circRNA}[0]} == 1); 
    for(my $i=0; $i < @{$circRNA{$circRNA}[0]}; $i++){
        if(exists $longest_tid{${$circRNA{$circRNA}[0]}[$i]}){
            my $length_tid_gid = $longest_tid{${$circRNA{$circRNA}[0]}[$i]}."_".${$circRNA{$circRNA}[0]}[$i];
#print $length_tid_gid,"\n";exit;
            push @{$circRNA{$circRNA}[1]},$length_tid_gid; 
            }    
        }
    @{$circRNA{$circRNA}[1]} = map {$_ ->[0]} sort {$b->[1] <=> $a->[1]} map {[$_,$_=~/(\d+)\_.+\_.+/]} @{$circRNA{$circRNA}[1]}; #sort by length, decreasing = T
    my @temp = split(/\_/,${$circRNA{$circRNA}[1]}[0]);
#print join"\,",@temp,"\n";exit;
    $circRNA{$circRNA}[2] = $temp[1]; 
    $circRNA{$circRNA}[3] = $temp[2];
#key = circRNAID, value[0] = gene_id, @{value[1]} = multi length_tid; value[2]=picked_tid, value[3]=picked_gid
    push @filtered_gene_ids,$circRNA{$circRNA}[3];
}

my %filtered_longest_tid;
foreach my $filtered_gene_id (@filtered_gene_ids){
    my @temp = split(/\_/,$longest_tid{$filtered_gene_id});
    $filtered_longest_tid{$filtered_gene_id."_".$temp[1]} = "";
}

#######################get the exon coordinates for the longest tid of each gene
my %exons;
open GTF,"<","$ARGV[0]"or die ($!);
while(<GTF>){
    chomp;
    next if($_ =~/^\#/); #skip description part of gtf
    my @line = split(/\t/,$_);
    if($line[2] eq "exon"){
        $line[8] =~ /gene_id \"(?<geneID>ENSG.+)\"\; transcript_id \"(?<transcriptID>ENST.+)\"\; gene_type \"(.+)\"\; gene_status \"(.+)\"\; gene_name \"(.+)\"\; transcript_type \"(.+)\"\; transcript_status \"(.+)\"\; transcript_name \"(.+)\"\; exon_number .*/;
        if (exists $filtered_longest_tid{$+{geneID}."_".$+{transcriptID}}){ # end > start
            push @{$exons{$+{geneID}}[0]},$line[3]; #exon start
            push @{$exons{$+{geneID}}[1]},$line[4]; #exon end
        }
    }
} 
close GTF;

#foreach  my $gene_ID (@gene_ids) {
#    print join"\,",@{$exons{$gene_ID}[0]},"\n";
#}
#exit
#######################print out
my $noExon_count;
open GENE,"<","$ARGV[1]"or die ($!);
while(<GENE>){
    chomp;
    my @line = split(/\t/,$_);
    my @range_start;
    my @range_end;
#    my @range_start_end; #remove duplicate 
    my @gene_id = split(/\,/,$line[7]);
    my @gene_name = split(/\,/,$line[6]);
    my $tID = $circRNA{$line[3]}[2];
    my $gID = $circRNA{$line[3]}[3];
    my $picked_geneName;
    for(my $i=0; $i<@gene_id; $i++){
        $picked_geneName = $gene_name[$i] if ($gene_id[$i] eq $gID);
    }
    for(my $i=0; $i<@{$exons{$gID}[0]}; $i++){ # end > start, no overlapped exons
        if(${$exons{$gID}[0]}[$i] >= $line[1] && ${$exons{$gID}[1]}[$i] <= $line[2]){ #within
            push @range_start,${$exons{$gID}[0]}[$i];
            push @range_end,${$exons{$gID}[1]}[$i];
#            push @range_start_end,${$exons{$line[7]}[0]}[$i]."_"."${$exons{$line[7]}[1]}[$i]";
        }elsif(${$exons{$gID}[0]}[$i] >= $line[1] && ${$exons{$gID}[0]}[$i] < $line[2] && ${$exons{$gID}[1]}[$i] > $line[2]){ #only once, right side overlap
            push @range_start,${$exons{$gID}[0]}[$i];
            push @range_end,$line[2];
        }elsif(${$exons{$gID}[0]}[$i] < $line[1] && ${$exons{$gID}[1]}[$i] > $line[1] && ${$exons{$gID}[1]}[$i] <= $line[2]){ #only once, left side overlap
            push @range_start,$line[1];
            push @range_end,${$exons{$gID}[1]}[$i];
        }elsif(${$exons{$gID}[0]}[$i] <= $line[1] && ${$exons{$gID}[1]}[$i] >= $line[2]){ #only once, cover
            push @range_start,$line[1];
            push @range_end,$line[2];
        }
    }
    my $starts = join"\,",@range_start;
    my $ends = join"\,",@range_end;
    if ($starts && $ends){
        push @line, ($starts,$ends,$tID,$gID,$picked_geneName); #insert before parameter 2 
    }else {
        $noExon_count += 1;
        push @line, ("-","-",$tID,$gID,$picked_geneName); #insert before parameter 2
    }
    print OUT join"\t",@line,"\n";
}
print "no exon: $noExon_count\n"; #80845
close GENE;
close OUT;
