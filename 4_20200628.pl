#usr/bin/perl -w
use strict;
print << "EOF";
This program is to rand(int) for IRES, m6A, ribosome profiling, and peptide composition.
EOF
die "Usage: perl $0  <input> <output> \n" unless @ARGV==2;

open IN,"<","$ARGV[0]"or die ($!);
open OUT,">","$ARGV[1]"or die ($!);

my $count;
while (<IN>){
    chomp;
    my @line = split(/\t/,$_);
    my $length = $line[2] - $line[1]; #or +1
    $count += 1;
    my $flag = $count % 3;
    if($flag == 2){
        my @IRES_start = (int(rand($length-1)), int(rand($length-1)), int(rand($length-1)));
        my $IRES_starts = join"\,",@IRES_start;
        my @IRES_end;
        for(my $i=0;$i<3;$i++){
            $IRES_end[$i] = $IRES_start[$i]+int(rand($length-$IRES_start[$i]));
        }
        my $IRES_ends = join"\,",@IRES_end;
#my @IRES_score = (rand(70),rand(70),rand(70));
        my $IRES_scores = rand(70);

        my $m6As = "-";
        my $ribo = "M";
        my $aa_composition = rand(1);
        
        my @temp = splice @line,14,0,($IRES_starts,$IRES_ends,$IRES_scores,$m6As,$ribo,$aa_composition); #insert before parameter 2
        print OUT join"\t",@line,"\n";
    }
    elsif($flag == 1){
        my @IRES_start = (int(rand($length-1)), int(rand($length-1)));
        my $IRES_starts = join"\,",@IRES_start;
        my @IRES_end;
        for(my $i=0;$i<2;$i++){
            $IRES_end[$i] = $IRES_start[$i]+int(rand($length-$IRES_start[$i]));
        }
        my $IRES_ends = join"\,",@IRES_end;
#my @IRES_score = (rand(50),rand(50));
        my $IRES_scores = rand(50);

        my @m6A = (int(rand($length-1)), int(rand($length-1)));
        my $m6As = join"\,",@m6A;
        my $ribo = "H";
        my $aa_composition = rand(1);
        
        my @temp = splice @line,14,0,($IRES_starts,$IRES_ends,$IRES_scores,$m6As,$ribo,$aa_composition); #insert before parameter 2
        print OUT join"\t",@line,"\n";
    }
    elsif($flag == 0){
        my $IRES_starts = "-";
        my $IRES_ends = "-";
        my $IRES_scores = 0;
        my @m6A = (int(rand($length-1)), int(rand($length-1)), int(rand($length-1)));
        my $m6As = join"\,",@m6A;
        my $ribo = "L";
        my $aa_composition = rand(1);

        my @temp = splice @line,14,0,($IRES_starts,$IRES_ends,$IRES_scores,$m6As,$ribo,$aa_composition); #insert before parameter 2
        print OUT join"\t",@line,"\n";
    }
}
close IN;
close OUT;
