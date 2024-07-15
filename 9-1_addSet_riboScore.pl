#!usr/bin/perl -w
use strict;
print << "EOF";
Add one independent ribosome/polysome-associated circRNA(BSJ) and re-calculate the score.
hg19 version
update the number of datasets each time

EOF

die "Usage: perl $0 <bsj.hg19> <PMID> <hg38_hg19_v2.0.tsv> <circ.RFP> <circ.newRFP>\n" unless @ARGV==5;

open BSJ, "<", "$ARGV[0]" or die ($!);
open LINK, "<", "$ARGV[2]" or die($!);
open RFP, "<", "$ARGV[3]" or die($!);
open OUT, ">", "$ARGV[4]" or die ($!);
########################################
my %BSJ;
while(<BSJ>){
    chomp;
    $BSJ{$_} = 0; #key = bsj.hg19
}
close BSJ;
#########################################
my %link;
<LINK>;
while(<LINK>){
    chomp;
    my @line = split(/\t/,$_);
    $link{$line[2]} = $line[3] if (exists $BSJ{$line[3]}); #key = bsj.hg38, value = bsj.hg19
}
close LINK;
#########################################
my $n = 3;

while(<RFP>){
    chomp;
    my @line = split(/\t/,$_);
    my @PMID = split(/\,/,$line[7]);
    if(exists $link{"$line[0]:$line[1]|$line[2]"} && !(grep {$ARGV[1] eq $_ } @PMID)){ 
        if ($PMID[0] ne "-"){
            push @PMID, $ARGV[1];
            $line[6] = @PMID/$n;
            $line[7] = join"\,",@PMID;
        }else{
            $line[6] = 1/$n;
            $line[7] = $ARGV[1];
        }
    }
    print OUT join"\t",@line,"\n";
}
close RFP;
close OUT;
