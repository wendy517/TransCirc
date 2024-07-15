#!usr/bin/perl -w
use strict;

print << "EOF";
This program is to calculate the frequency of di-aa (amino acid) within one peptide.
EOF

die "Usage: perl $0 <in.pep> <out.pep.freq-diaa>  \n" unless @ARGV==2;

open IN, "<", "$ARGV[0]" or die ($!);
open OUT, ">", "$ARGV[1]" or die ($!);

my $aa = "ARNDCQEGHILKMFPSTWYV";
my @aa = split(//,$aa);
my %di_aa;
for(my $i=0; $i<@aa; $i++){
    for(my $j=0; $j<@aa; $j++){
        $di_aa{$aa[$i].$aa[$j]} = 0;        
    }
}

my @di_aa = sort keys %di_aa;
my @count = values %di_aa;
my $di_aa = join"\t",@di_aa;
my $count = join"\t",@count;
print OUT "chrom\tcircStart\tcircEnd\tcircID\t.\tstrand\tframe\torfStart\torfEnd\tacross\torfSeq\tpepSeq\t$di_aa\n";

while(<IN>){
    chomp;
    next if ($_ =~ /^\#/);
    my @line = split(/\t/,$_);
    print OUT $_.$count,"\n" if ($line[-1] eq "-");
    next if ($line[-1] eq "-");

    foreach my $diaa (@di_aa){$di_aa{$diaa} = 0;} #initiation

    my $bare_pep = substr($line[-1], 1, length($line[-1])-2); ##ignor the start and the stop codon
    for(my $i=0; $i<length($bare_pep)-1; $i++){ 
        my $di_aa = substr($bare_pep, $i, 2);
        $di_aa{$di_aa} += 1;
    }   
    
    my @freq;
    foreach my $diaa (@di_aa){
        push @freq, ($di_aa{$diaa}/(length($bare_pep)-1));
    }
    my $freq_di_aa = join"\t",@freq;
    print OUT $_.$freq_di_aa,"\n";
}
close IN;
close OUT;
