#!usr/bin/perl -w
use strict;

print << "EOF";
This program is to extract potential ORFs in circRNA spliced sequences and provide translated peptide sequences accordingly.
Principals:
0. All potential ORFs must across BSJ and have minimal length of 20aa(63nt for cORF, 60nt for rcORF).
1. Alternative start codons[NTG]: ATG > CTG > GTG > TTG
2. Pick the longest ORF within one frame.
3. Priority:1>2
*Make sure input sequences are in stranded order.
EOF
die "Usage: perl $0 <in.seq> <out.ORF.row> <out.ORF.column> <badCondon.txt> \n" unless @ARGV==4;

open IN, "<", "$ARGV[0]" or die ($!);
open OUTR, ">", "$ARGV[1]" or die ($!);
open OUTC, ">", "$ARGV[2]" or die ($!);
open BAD, ">", "$ARGV[3]" or die ($!);

my @alternative_start_codons = qw(ATG CTG GTG TTG);
my @stop_codons = qw(TAG TAA TGA);
my %hash_stop_codons = map{$_=>1} @stop_codons;

######statistics
my $total_input = 0;
my $total_seq = 0;
my $ORF = 0; #0 < x < 3*$total_seq
my $no_ORF = 0; #x+$ORF = 3*$total_seq
my $cORF = 0; # 0< x < $ORF
my $rcORF = 0; # 0< x < $rcORF, x+$cORF = $ORF
my $triORF_seq = 0; #0< x < $total_seq
my $diORF_seq = 0; #0< x < $total_seq
my $singleORF_seq = 0; #0< x < $total_seq
my $noneORF_seq = 0; #0< x < $total_seq, x+ $singleORF_seq + $diORF_seq + $triORF_seq = $total_seq
my $rcORFcontain_seq = 0;
my $no_start_ORF = 0; #0< x < 3*$noneORF_seq
my $less_than20aa_ORF = 0; #0< x < 3*$noneORF_seq
my $across0_ORF = 0; #0< x < 3*$noneORF_seq, x+ $no_start_ORF  + $less_than20aa_ORF = $no_ORF 
my $badPeptide = 0;

#######main
while(<IN>){
    chomp;
    $total_input += 1;
    my $tag = 0;
    my @line = split(/\t/,$_);
    my @circ = @line[0..5];
    if ($line [13] eq "-"){
        print OUTR join"\t",(@line,"-","-","-","-","-", "-","-","-","-","-", "-","-","-","-","-"),"\n";
        print OUTC join"\t",(@circ,"frame0","-","-","-","-","-"),"\n";
        print OUTC join"\t",(@circ,"frame1","-","-","-","-","-"),"\n";
        print OUTC join"\t",(@circ,"frame2","-","-","-","-","-"),"\n";
    }
    next if ($line[13] eq "-");
    $total_seq += 1;

    my $frame0_start; my $frame0_stop; my $frame0_across; my $frame0_ORF; my $frame0_pep; my $frame0_noORF;
    my $frame1_start; my $frame1_stop; my $frame1_across; my $frame1_ORF; my $frame1_pep; my $frame1_noORF;
    my $frame2_start; my $frame2_stop; my $frame2_across; my $frame2_ORF; my $frame2_pep; my $frame2_noORF;
    my $frame0_noORF = "-"; my $frame1_noORF = "-"; my $frame2_noORF = "-";
    
    ### seq in-frame
    if (length($line[13])%3 == 0){ 
        my @nt = split(//,$line[13]); #examine the 2 codons at BSJ
        my $seq_BSJ = $line[13].$nt[0].$nt[1];
        my @frame0_codons; my @frame1_codons; my @frame2_codons;
        for(my $i=0; $i<=length($seq_BSJ)-3; $i+=3){
            my $codon = substr($seq_BSJ,$i,3);
            push @frame0_codons,$codon;
        }
        for(my $i=1; $i<=length($seq_BSJ)-3; $i+=3){
            my $codon = substr($seq_BSJ,$i,3);
            push @frame1_codons,$codon;
        }
        for(my $i=2; $i<=length($seq_BSJ)-3; $i+=3){
            my $codon = substr($seq_BSJ,$i,3);
            push @frame2_codons,$codon;
        }
        ###in-frame:frame0
        ($frame0_start, $frame0_stop, $frame0_across, $frame0_ORF, $frame0_noORF) = inFrame_ORF(\@frame0_codons, 0, $line[13], \@stop_codons);
        $frame0_pep = translation($frame0_ORF);
        push @line, ($frame0_start,$frame0_stop,$frame0_across,$frame0_ORF,$frame0_pep);
        my @circ_frame0=(@circ,"frame0",$frame0_start,$frame0_stop,$frame0_across,$frame0_ORF,$frame0_pep);
        print OUTC join"\t",@circ_frame0,"\n";
        ###in-frame:frame1
        ($frame1_start, $frame1_stop, $frame1_across, $frame1_ORF, $frame1_noORF) = inFrame_ORF(\@frame1_codons, 1, $line[13], \@stop_codons);
        $frame1_pep = translation($frame1_ORF);
        push @line, ($frame1_start,$frame1_stop,$frame1_across,$frame1_ORF,$frame1_pep);
        my @circ_frame1=(@circ,"frame1",$frame1_start,$frame1_stop,$frame1_across,$frame1_ORF,$frame1_pep);
        print OUTC join"\t",@circ_frame1,"\n";
        ###in-frame:frame2
        ($frame2_start, $frame2_stop, $frame2_across, $frame2_ORF, $frame2_noORF) = inFrame_ORF(\@frame2_codons, 2, $line[13], \@stop_codons);
        $frame2_pep = translation($frame2_ORF);
        push @line, ($frame2_start,$frame2_stop,$frame2_across,$frame2_ORF,$frame2_pep);
        my @circ_frame2=(@circ,"frame2",$frame2_start,$frame2_stop,$frame2_across,$frame2_ORF,$frame2_pep);
        print OUTC join"\t",@circ_frame2,"\n";
    }

    ###seq not in-frame
    else{ 
        my $qua_seq = $line[13] x 4;
        my @qua_seq_codons;
        for (my $i = 0; $i <=length($qua_seq)-3; $i+=3){
            my $codon = substr($qua_seq, $i, 3);
            push @qua_seq_codons,$codon;
        }
#print join"\,",@qua_seq_codons,"\n";exit;
        ###initialization
        my @frame0_startCodon= (); my @frame1_startCodon= (); my @frame2_startCodon= () ;# two-dimensional array,c1=4, c2=index of each start codon in codons array
        @{$frame0_startCodon[0]}= (); @{$frame1_startCodon[0]}= (); @{$frame2_startCodon[0]}= ();
        @{$frame0_startCodon[1]}= (); @{$frame1_startCodon[1]}= (); @{$frame2_startCodon[1]}= ();
        @{$frame0_startCodon[2]}= (); @{$frame1_startCodon[2]}= (); @{$frame2_startCodon[2]}= ();
        @{$frame0_startCodon[3]}= (); @{$frame1_startCodon[3]}= (); @{$frame2_startCodon[3]}= ();
        my @quaSeq_stopCodon = ();
        ###codon block 0
        for(my $i = 0; $i <= (length($line[13])-1)/3 ; $i ++){
            if ($qua_seq_codons[$i] eq "ATG"){push @{$frame0_startCodon[0]},$i;}
            elsif($qua_seq_codons[$i] eq "CTG"){push @{$frame0_startCodon[1]},$i;}
            elsif($qua_seq_codons[$i] eq "GTG"){push @{$frame0_startCodon[2]},$i;}
            elsif($qua_seq_codons[$i] eq "TTG"){push @{$frame0_startCodon[3]},$i;}        
            elsif($qua_seq_codons[$i] eq "TAG" || $qua_seq_codons[$i] eq "TAA" || $qua_seq_codons[$i] eq "TGA"){push @quaSeq_stopCodon, $i; }
        }
        ###codon block 1
        for (my $i = int((length($line[13])-1)/3) + 1 ; $i <= (2*length($line[13]) -1)/3; $i ++){
            if(length($line[13]) % 3 == 1){
                if ($qua_seq_codons[$i] eq "ATG"){push @{$frame2_startCodon[0]},$i;}
                elsif($qua_seq_codons[$i] eq "CTG"){push @{$frame2_startCodon[1]},$i;}
                elsif($qua_seq_codons[$i] eq "GTG"){push @{$frame2_startCodon[2]},$i;}
                elsif($qua_seq_codons[$i] eq "TTG"){push @{$frame2_startCodon[3]},$i;}
                elsif($qua_seq_codons[$i] eq "TAG" || $qua_seq_codons[$i] eq "TAA" || $qua_seq_codons[$i] eq "TGA"){push @quaSeq_stopCodon, $i; }
            }
            elsif (length($line[13]) % 3 == 2){
                if ($qua_seq_codons[$i] eq "ATG"){push @{$frame1_startCodon[0]},$i;}
                elsif($qua_seq_codons[$i] eq "CTG"){push @{$frame1_startCodon[1]},$i;}
                elsif($qua_seq_codons[$i] eq "GTG"){push @{$frame1_startCodon[2]},$i;}
                elsif($qua_seq_codons[$i] eq "TTG"){push @{$frame1_startCodon[3]},$i;}
                elsif($qua_seq_codons[$i] eq "TAG" || $qua_seq_codons[$i] eq "TAA" || $qua_seq_codons[$i] eq "TGA"){push @quaSeq_stopCodon, $i; }
            }
        }
        ###codon block 2
        for (my $i = int((2*length($line[13])-1)/3) + 1; $i <= (3*length($line[13]) -1)/3; $i ++){
            if(length($line[13]) % 3 == 1){
                if ($qua_seq_codons[$i] eq "ATG"){push @{$frame1_startCodon[0]},$i;}
                elsif($qua_seq_codons[$i] eq "CTG"){push @{$frame1_startCodon[1]},$i;}
                elsif($qua_seq_codons[$i] eq "GTG"){push @{$frame1_startCodon[2]},$i;}
                elsif($qua_seq_codons[$i] eq "TTG"){push @{$frame1_startCodon[3]},$i;}
                elsif($qua_seq_codons[$i] eq "TAG" || $qua_seq_codons[$i] eq "TAA" || $qua_seq_codons[$i] eq "TGA"){push @quaSeq_stopCodon, $i; }
            }
            elsif (length($line[13]) % 3 == 2){
                if ($qua_seq_codons[$i] eq "ATG"){push @{$frame2_startCodon[0]},$i;}
                elsif($qua_seq_codons[$i] eq "CTG"){push @{$frame2_startCodon[1]},$i;}
                elsif($qua_seq_codons[$i] eq "GTG"){push @{$frame2_startCodon[2]},$i;}
                elsif($qua_seq_codons[$i] eq "TTG"){push @{$frame2_startCodon[3]},$i;}
                elsif($qua_seq_codons[$i] eq "TAG" || $qua_seq_codons[$i] eq "TAA" || $qua_seq_codons[$i] eq "TGA"){push @quaSeq_stopCodon, $i; }
            }
        }
        ###codon block 3
        for (my $i = length($line[13]); $i < @qua_seq_codons; $i ++){
            if($qua_seq_codons[$i] eq "TAG" || $qua_seq_codons[$i] eq "TAA" || $qua_seq_codons[$i] eq "TGA"){push @quaSeq_stopCodon, $i; }
        }

#print join "\t",@{$frame1_startCodon[2]},"\n";exit;
        ###non-in-frame seq: examine rcORF for each frame
        if (!(@quaSeq_stopCodon)){ 
            ($frame0_start, $frame0_stop, $frame0_across, $frame0_ORF, $frame0_noORF) = non_inFrame_rcORF(\@frame0_startCodon, $line[13]);
            $frame0_pep = translation($frame0_ORF);
            push @line, ($frame0_start, $frame0_stop, $frame0_across, $frame0_ORF, $frame0_pep);
            my @circ_frame0=(@circ,"frame0", $frame0_start, $frame0_stop, $frame0_across, $frame0_ORF, $frame0_pep);
            print OUTC join"\t",@circ_frame0,"\n";
            
            ($frame1_start, $frame1_stop, $frame1_across, $frame1_ORF, $frame1_noORF) = non_inFrame_rcORF(\@frame1_startCodon, $line[13]);
            $frame1_pep = translation($frame1_ORF);
            push @line, ($frame1_start, $frame1_stop, $frame1_across, $frame1_ORF, $frame1_pep);
            my @circ_frame1=(@circ,"frame1", $frame1_start, $frame1_stop, $frame1_across, $frame1_ORF, $frame1_pep);
            print OUTC join"\t",@circ_frame1,"\n";

            ($frame2_start, $frame2_stop, $frame2_across, $frame2_ORF, $frame2_noORF) = non_inFrame_rcORF(\@frame2_startCodon, $line[13]);
            $frame2_pep = translation($frame2_ORF);
            push @line, ($frame2_start, $frame2_stop, $frame2_across, $frame2_ORF, $frame2_pep);
            my @circ_frame2=(@circ,"frame2", $frame2_start, $frame2_stop, $frame2_across, $frame2_ORF, $frame2_pep);
            print OUTC join"\t",@circ_frame2,"\n";
        }
        ###non-in-frame seq: examine cORF for each frame
        else{ 
            ($frame0_start, $frame0_stop, $frame0_across, $frame0_ORF, $frame0_noORF) = non_inFrame_cORF(\@frame0_startCodon, \@quaSeq_stopCodon, $line[13]);
            $frame0_pep = translation($frame0_ORF);
            push @line, ($frame0_start, $frame0_stop, $frame0_across, $frame0_ORF, $frame0_pep);
            my @circ_frame0=(@circ,"frame0", $frame0_start, $frame0_stop, $frame0_across, $frame0_ORF, $frame0_pep);
            print OUTC join"\t",@circ_frame0,"\n";

            ($frame1_start, $frame1_stop, $frame1_across, $frame1_ORF, $frame1_noORF) = non_inFrame_cORF(\@frame1_startCodon, \@quaSeq_stopCodon, $line[13]);
            $frame1_pep = translation($frame1_ORF);
            push @line, ($frame1_start, $frame1_stop, $frame1_across, $frame1_ORF, $frame1_pep);
            my @circ_frame1=(@circ,"frame1", $frame1_start, $frame1_stop, $frame1_across, $frame1_ORF, $frame1_pep);
            print OUTC join"\t",@circ_frame1,"\n";

            ($frame2_start, $frame2_stop, $frame2_across, $frame2_ORF, $frame2_noORF) = non_inFrame_cORF(\@frame2_startCodon, \@quaSeq_stopCodon, $line[13]);
            $frame2_pep = translation($frame2_ORF);
            push @line, ($frame2_start, $frame2_stop, $frame2_across, $frame2_ORF, $frame2_pep);
            my @circ_frame2=(@circ,"frame2", $frame2_start, $frame2_stop, $frame2_across, $frame2_ORF, $frame2_pep);
            print OUTC join"\t",@circ_frame2,"\n";
        }
    }

#######print & statistic
    print OUTR join"\t",@line,"\n";
    my @across = ($frame0_across, $frame1_across, $frame2_across);
    $rcORFcontain_seq += 1 if grep {"3" eq $_} @across;    
    foreach my $frame_across (@across){
        $ORF += 1 if ($frame_across ne "-");
        $tag += 1 if ($frame_across ne "-");
        $no_ORF += 1 if ($frame_across eq "-");
        next if ($frame_across eq "-");
        $cORF += 1 if ($frame_across == 1 || $frame_across == 2);
        $rcORF += 1 if ($frame_across == 3);
    }
    my @noORF_type = ($frame0_noORF, $frame1_noORF, $frame2_noORF);
    foreach my $frame_noORF (@noORF_type){
        $no_start_ORF += 1 if ($frame_noORF eq "no_start");
        $less_than20aa_ORF += 1 if ($frame_noORF eq "less_than20aa");
        $across0_ORF += 1 if ($frame_noORF eq "across0");
    }
    my @peptide = ($frame0_pep, $frame1_pep, $frame2_pep);
    foreach my $frame_pep (@peptide){
        $badPeptide += 1 if ($frame_pep =~ /X/);
    }
    $triORF_seq += 1 if ($tag == 3);
    $diORF_seq += 1 if ($tag == 2);
    $singleORF_seq += 1 if ($tag == 1);
    $noneORF_seq += 1 if ($tag == 0);
}
close IN;
close OUTR;
close OUTC;
close BAD;

print "$total_input total input\n";
print "$total_seq input with sequence\n";
print "$ORF total ORFs\n";
print "$no_ORF failed ORFs\n";
print "$cORF total cORFs\n";
print "$rcORF total rcORFs\n";
print "$triORF_seq sequences with 3 ORFs\n";
print "$diORF_seq sequences with 2 ORFs\n";
print "$singleORF_seq sequences with 1 ORFs\n";
print "$noneORF_seq sequences with no ORFs\n";
print "$rcORFcontain_seq sequences contain rcORFs\n";
print "$no_start_ORF failed ORFs without start codon(s)\n";
print "$across0_ORF failed ORFs across BSJ 0 times\n";
print "$less_than20aa_ORF failed ORFs length < 60nt (20aa)\n";
print "$badPeptide peptides have uncertain codon\n";


###################################################################################
sub inFrame_ORF{ 
#input: (\@frame_codons, $frame, $line[13], \@stop_codons)
#output: ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF)
    my ($ref_frame_codons, $frame, $seq, $ref_stop_codons) = @_;
    my @frame_codons = @$ref_frame_codons;
    my @stop_codons = @$ref_stop_codons;

    my @frame_startCodon;
    @{$frame_startCodon[0]}= (); 
    @{$frame_startCodon[1]}= ();
    @{$frame_startCodon[2]}= ();
    @{$frame_startCodon[3]}= ();
    my $frame_start; my $frame_stop; my $frame_across; my $frame_ORF; my $frame_noORF = "-";

    for(my $i=0; $i<@frame_codons; $i++){
        if ($frame_codons[$i] eq "ATG"){ push @{$frame_startCodon[0]}, $i;}
        elsif ($frame_codons[$i] eq "CTG"){ push @{$frame_startCodon[1]}, $i;}
        elsif ($frame_codons[$i] eq "GTG"){ push @{$frame_startCodon[2]}, $i;}
        elsif ($frame_codons[$i] eq "TTG"){ push @{$frame_startCodon[3]}, $i;}
    }

    my %hash_stop_codons = map{$_=>1} @stop_codons;
    my @common_stop = grep {$hash_stop_codons{$_}} @frame_codons;
    ###examine no stop (across = 3 ?)
    if (!(@common_stop)){ 
        if (!(@{$frame_startCodon[0]} || @{$frame_startCodon[1]} || @{$frame_startCodon[2]} || @{$frame_startCodon[3]})){
            $frame_noORF = "no_start";
            $frame_start = "-"; $frame_stop = "-"; $frame_across = "-"; $frame_ORF = "-";
            return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
        }
        if (@{$frame_startCodon[0]}){ #pick any start codon(the first one)
            $frame_start = 3 * ${$frame_startCodon[0]}[0] + $frame;
            $frame_stop = $frame_start - 1;
            $frame_across = 3;
            my $frame_str_a = substr($seq, $frame_start);#from the start to the right BSJ
            my $frame_str_b = substr($seq, 0, $frame_stop + 1); #from the left BSJ to the end
            $frame_ORF = ($frame_str_a.$frame_str_b) x $frame_across; #repeat 3 times
            if(length($frame_ORF) >= 60){
                $frame_noORF = "-";
                return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
            }else{
                $frame_ORF = "";
                $frame_noORF = "less_than20aa";
            }
        } #examine length quality first, then alternative start codon
        if (!($frame_ORF) && @{$frame_startCodon[1]}){
            $frame_start = 3 * ${$frame_startCodon[1]}[0] + $frame;
            $frame_stop = $frame_start - 1;
            $frame_across = 3;
            my $frame_str_a = substr($seq, $frame_start);
            my $frame_str_b = substr($seq, 0, $frame_stop + 1);
            $frame_ORF = ($frame_str_a.$frame_str_b) x $frame_across;
            if(length($frame_ORF) >= 60){
                $frame_noORF = "-";
                return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
            }else{
                $frame_ORF = "";
                $frame_noORF = "less_than20aa";
            }
        }
        if (!($frame_ORF) && @{$frame_startCodon[2]}){
            $frame_start = 3 * ${$frame_startCodon[2]}[0] + $frame;
            $frame_stop = $frame_start - 1;
            $frame_across = 3;
            my $frame_str_a = substr($seq, $frame_start);
            my $frame_str_b = substr($seq, 0, $frame_stop + 1);
            $frame_ORF = ($frame_str_a.$frame_str_b) x $frame_across;
            if(length($frame_ORF) >= 60){
                $frame_noORF = "-";
                return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
            }else{
                $frame_ORF = "";
                $frame_noORF = "less_than20aa";
            }
        }
        if (!($frame_ORF) && @{$frame_startCodon[3]}){
            $frame_start = 3 * ${$frame_startCodon[3]}[0] + $frame;
            $frame_stop = $frame_start - 1;
            $frame_across = 3;
            my $frame_str_a = substr($seq, $frame_start);
            my $frame_str_b = substr($seq, 0, $frame_stop + 1);
            $frame_ORF = ($frame_str_a.$frame_str_b) x $frame_across;
            if(length($frame_ORF) >= 60){
                $frame_noORF = "-";
                return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
            }else{
                $frame_ORF = "";
                $frame_noORF = "less_than20aa";
            }
        }
        if (!($frame_ORF) && @{$frame_startCodon[0]} || @{$frame_startCodon[1]} || @{$frame_startCodon[2]} || @{$frame_startCodon[3]}){
            return ("-","-","-","-", $frame_noORF);
        }
    }
    ###examine have stop codon(across = 0 or 1 ?)
    else{
        #examine start codon first
        if (!(@{$frame_startCodon[0]} || @{$frame_startCodon[1]} || @{$frame_startCodon[2]} || @{$frame_startCodon[3]})){
            $frame_noORF = "no_start";
            $frame_start = "-"; $frame_stop = "-"; $frame_across = "-"; $frame_ORF = "-";
            return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
        } 
        else{ #have start codon (and stop codon)
            if ($frame == 0){
                $frame_noORF = "across0";
                $frame_start = "-"; $frame_stop = "-"; $frame_across = "-"; $frame_ORF = "-";
                return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
            }
            if ($frame == 1 || $frame == 2){
                my @rest_condons = @frame_codons[0..$#frame_codons-1]; 
                my @rest_common_stop = grep {$hash_stop_codons{$_}} @rest_condons;
                ###only if the last condon(at BSJ) is stop codon   
                if(@rest_common_stop){
                    $frame_noORF = "across0";
                    $frame_start = "-"; $frame_stop = "-"; $frame_across = "-"; $frame_ORF = "-";
                    return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
                }
                elsif ($frame_codons[-1] ~~ @stop_codons) {
                    $frame_stop = $frame - 1; #stop_right
                    if (@{$frame_startCodon[0]}){
                        $frame_start = 3 * ${$frame_startCodon[0]}[0] + $frame;
                        $frame_across = 1;
                        my $frame_str_a = substr($seq, $frame_start);
                        my $frame_str_b = substr($seq, 0, $frame_stop + 1);
                        $frame_ORF = ($frame_str_a.$frame_str_b);
                        if(length($frame_ORF) >= 63){
                            $frame_noORF = "-";
                            return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
                        }else{
                            $frame_ORF = "";
                            $frame_noORF = "less_than20aa";
                        }
                    }
                    if (!($frame_ORF) && @{$frame_startCodon[1]}){
                        $frame_start = 3 * ${$frame_startCodon[1]}[0] + $frame;
                        $frame_across = 1;
                        my $frame_str_a = substr($seq, $frame_start);
                        my $frame_str_b = substr($seq, 0, $frame_stop + 1);
                        $frame_ORF = ($frame_str_a.$frame_str_b);
                        if(length($frame_ORF) >= 63){
                            $frame_noORF = "-";
                            return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
                        }else{
                            $frame_ORF = "";
                            $frame_noORF = "less_than20aa";
                        }
                    }
                    if (!($frame_ORF) && @{$frame_startCodon[2]}){
                        $frame_start = 3 * ${$frame_startCodon[2]}[0] + $frame;
                        $frame_across = 1;
                        my $frame_str_a = substr($seq, $frame_start);
                        my $frame_str_b = substr($seq, 0, $frame_stop + 1);
                        $frame_ORF = ($frame_str_a.$frame_str_b);
                        if(length($frame_ORF) >= 63){
                            $frame_noORF = "-";
                            return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
                        }else{
                            $frame_ORF = "";
                            $frame_noORF = "less_than20aa";
                        }
                    }
                    if (!($frame_ORF) && @{$frame_startCodon[3]}){
                        $frame_start = 3 * ${$frame_startCodon[3]}[0] + $frame;
                        $frame_across = 1;
                        my $frame_str_a = substr($seq, $frame_start);
                        my $frame_str_b = substr($seq, 0, $frame_stop + 1);
                        $frame_ORF = ($frame_str_a.$frame_str_b);
                        if(length($frame_ORF) >= 63){
                            $frame_noORF = "-";
                            return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
                        }else{
                            $frame_ORF = "";
                            $frame_noORF = "less_than20aa";
                        }
                    }
                    if (!($frame_ORF)){ #
                        return ("-","-","-","-", $frame_noORF);
                    }
                }
            }
        }
    }
}

###########################################################################################
sub non_inFrame_rcORF{
#input:(\@frame_startCodon, $seq)
#output: ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF)
    my ($ref, $seq) = @_;
    my @frame_startCodon = @$ref;
    my $frame_start; my $frame_stop; my $frame_across; my $frame_ORF; my $frame_noORF = "-";

    if (!(@{$frame_startCodon[0]} || @{$frame_startCodon[1]} || @{$frame_startCodon[2]} || @{$frame_startCodon[3]})){
        $frame_noORF = "no_start";
        $frame_start = "-"; $frame_stop = "-"; $frame_across = "-"; $frame_ORF = "-"; 
        return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
    }
    if(@{$frame_startCodon[0]}){
        my $start = 3 * ${$frame_startCodon[0]}[0];
        if ($start < length($seq)) {
            $frame_start = $start;
        }elsif($start < 2*length($seq)){
            $frame_start = $start - length($seq);
        }elsif($start < 3*length($seq)){
            $frame_start = $start - 2*length($seq);
        }else{
            $frame_start = $start - 3*length($seq);
        }
        $frame_stop = $frame_start - 1;
        $frame_across = 3;
        my $frame_str_a = substr($seq, $frame_start); 
        my $frame_str_b = substr($seq, 0, $frame_stop + 1);
        $frame_ORF = ($frame_str_a.$frame_str_b) x $frame_across;
        if(length($frame_ORF) >= 60){
            $frame_noORF = "-";
            return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
        }else {
            $frame_ORF = ""; 
            $frame_noORF = "less_than20aa";
        }
    } #examine length quality first, then alternative start codon
    if(!($frame_ORF) && @{$frame_startCodon[1]}){ 
        my $start = 3 * ${$frame_startCodon[1]}[0];
        if ($start < length($seq)){
            $frame_start = $start;
        }elsif($start < 2*length($seq)){
            $frame_start = $start - length($seq);
        }elsif($start < 3*length($seq)){
            $frame_start = $start - 2*length($seq);
        }else{
            $frame_start = $start - 3*length($seq);
        }
        $frame_stop = $frame_start - 1;
        $frame_across = 3;
        my $frame_str_a = substr($seq, $frame_start);
        my $frame_str_b = substr($seq, 0, $frame_stop + 1);
        $frame_ORF = ($frame_str_a.$frame_str_b) x $frame_across;
        if(length($frame_ORF) >= 60){
            $frame_noORF = "-";
            return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
        }else {
            $frame_ORF = "";
            $frame_noORF = "less_than20aa";
        }
    }
    if(!($frame_ORF) && @{$frame_startCodon[2]}){
        my $start = 3 * ${$frame_startCodon[2]}[0];
        if ($start < length($seq)){
            $frame_start = $start;
        }elsif($start < 2*length($seq)){
            $frame_start = $start - length($seq);
        }elsif($start < 3*length($seq)){
            $frame_start = $start - 2*length($seq);
        }else{
            $frame_start = $start - 3*length($seq);
        }
        $frame_stop = $frame_start - 1;
        $frame_across = 3;
        my $frame_str_a = substr($seq, $frame_start);
        my $frame_str_b = substr($seq, 0, $frame_stop + 1);
        $frame_ORF = ($frame_str_a.$frame_str_b) x $frame_across;
        if(length($frame_ORF) >= 60){
            $frame_noORF = "-";
            return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
        }else{
            $frame_ORF = "";
            $frame_noORF = "less_than20aa";
        }
    }
    if(!($frame_ORF) && @{$frame_startCodon[3]}){
        my $start = 3 * ${$frame_startCodon[3]}[0];
        if ($start < length($seq)){
            $frame_start = $start;
        }elsif($start < 2*length($seq)){
            $frame_start = $start - length($seq);
        }elsif($start < 3*length($seq)){
            $frame_start = $start - 2*length($seq);
        }else{
            $frame_start = $start - 3*length($seq);
        }
        $frame_stop = $frame_start - 1;
        $frame_across = 3;
        my $frame_str_a = substr($seq, $frame_start);
        my $frame_str_b = substr($seq, 0, $frame_stop + 1);
        $frame_ORF = ($frame_str_a.$frame_str_b) x $frame_across;
        if(length($frame_ORF) >= 60){
            $frame_noORF = "-";
            return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
        }else{
            $frame_ORF = "";
            $frame_noORF = "less_than20aa";
        }
    }
    
   if (!($frame_ORF) && (@{$frame_startCodon[0]} || @{$frame_startCodon[1]} || @{$frame_startCodon[2]} || @{$frame_startCodon[3]})){ 
        return ("-","-","-","-", $frame_noORF);
    }    
}

###########################################################################################
sub non_inFrame_cORF{
#input:(\@frame_startCodon, \@quaSeq_stopCodon, $seq)
#output: ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF)
    my ($ref, $ref_stopCodon, $seq) = @_;
    my @frame_startCodon = @$ref;
    my @quaSeq_stopCodon = @$ref_stopCodon;
    
    my $frame_start; my $frame_stop; my $frame_across; my $frame_ORF; my $frame_noORF = "-";
    
    my @frame_start_stop;
    @{$frame_start_stop[0]} = (); 
    @{$frame_start_stop[1]} = ();
    @{$frame_start_stop[2]} = (); 
    @{$frame_start_stop[3]} = ();
    
    if (!(@{$frame_startCodon[0]} || @{$frame_startCodon[1]} || @{$frame_startCodon[2]} || @{$frame_startCodon[3]})){ #no start codon in the frame
        $frame_noORF = "no_start";
        $frame_start = "-"; $frame_stop = "-"; $frame_across = "-"; $frame_ORF = "-";
        return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
    }


    if(@{$frame_startCodon[0]}){
        my %hash_start;
        START:foreach my $start (@{$frame_startCodon[0]}){
                foreach my $stop (@quaSeq_stopCodon){
                    ### the stop codon after the start codon can not within one block
                    next if ($stop < $start);
                    if (3 * $start < length($seq)){
                        next START if ((3*$stop+2 < length($seq)) && (3*$stop+2 > 3*$start));
                        my $distance = (3*$stop+2) - (3*$start);
                        push @{$hash_start{$start}}, $distance."_".$start."_".$stop;
                    }
                    elsif (3 * $start < 2 * length($seq)){
                        next START if ((3*$stop+2 < 2 * length($seq)) && (3*$stop+2 > 3*$start));
                        my $distance = (3*$stop+2) - (3*$start);
                        push @{$hash_start{$start}}, $distance."_".$start."_".$stop;
                    }
                    elsif (3 * $start < 3 * length($seq)){
                        next START if ((3*$stop+2 < 3 * length($seq)) && (3*$stop+2 > 3*$start));
                        my $distance = (3*$stop+2) - (3*$start);
                        push @{$hash_start{$start}}, $distance."_".$start."_".$stop;
                    } # no start codon in block 3
                }
        }
        foreach my $start (@{$frame_startCodon[0]}){
            next unless (exists $hash_start{$start});
            @{$hash_start{$start}} = map {$_ ->[0]} sort {$a->[1] <=> $b->[1]} map {[$_,$_=~/(\d+)\_.+\_.+/]} @{$hash_start{$start}}; #sort by distance, decreasing = F, pick the cloest stop codon for each start codon
            push @{$frame_start_stop[0]}, ${$hash_start{$start}}[0];
        }
        
    if (!@{$frame_start_stop[0]}){ $frame_noORF = "across0";}
        else{
            @{$frame_start_stop[0]} = map {$_ ->[0]} sort {$b->[1] <=> $a->[1]} map {[$_,$_=~/(\d+)\_.+\_.+/]} @{$frame_start_stop[0]}; #sort by distance, decreasing = T, pick the longest ORF for one startCdon in each frame
            my @temp = split(/\_/,${$frame_start_stop[0]}[0]);
            my $start = 3 * $temp[1]; my $stop = 3 * $temp[2] + 2;
#print $start,"\t",$stop,"\n";exit;
            if ($start < length($seq)){
                $frame_start = $start;
                if ($stop < 2*length($seq)){
                    $frame_stop = $stop - length($seq);
                    $frame_across = 1;
                }
                elsif ($stop < 3*length($seq)){
                    $frame_stop = $stop - 2*length($seq);
                    $frame_across = 2;
                }
                else { 
                   $frame_stop = $stop - 3*length($seq); 
                   $frame_across = 1;
                }
            }
            elsif($start < 2*length($seq)){
                $frame_start = $start - length($seq);
                if ($stop < 3*length($seq)){
                    $frame_stop = $stop - 2*length($seq);
                    $frame_across = 1;
                }
                else{
                    $frame_stop = $stop - 3*length($seq);
                    $frame_across = 2;
                }
            }
            elsif($start < 3*length($seq)){
                $frame_start = $start - 2*length($seq);
                $frame_stop = $stop - 3*length($seq);
                $frame_across = 1;
            }
#print $frame_start,"\t",$frame_stop,"\n";exit;            
            my $frame_str_a = substr($seq, $frame_start);
            my $frame_str_b = substr($seq, 0, $frame_stop + 1);
            if ($frame_across == 1){ $frame_ORF = $frame_str_a.$frame_str_b;}
            elsif ($frame_across == 2){$frame_ORF = $frame_str_a.$seq.$frame_str_b;}
#else{ $frame_ORF = "";}

            if(length($frame_ORF) >= 63){
                $frame_noORF = "-";
                return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
            }else{
                $frame_ORF = "";
                $frame_noORF = "less_than20aa";
            }
        }
    }
    
    
    if(!($frame_ORF) && @{$frame_startCodon[1]}){
#print join"\t",@{$frame_startCodon[1]},"\n";  
        my %hash_start;
        START:foreach my $start (@{$frame_startCodon[1]}){
                foreach my $stop (@quaSeq_stopCodon){
                    next if ($stop < $start);
                    if (3 * $start < length($seq)){
                        next START if ((3*$stop+2 < length($seq)) && (3*$stop+2 > 3*$start));
                        my $distance = (3*$stop+2) - (3*$start);
                        push @{$hash_start{$start}}, $distance."_".$start."_".$stop;
                    }
                    elsif (3 * $start < 2 * length($seq)){
                        next START if ((3*$stop+2 < 2 * length($seq)) && (3*$stop+2 > 3*$start));
                        my $distance = (3*$stop+2) - (3*$start);
                        push @{$hash_start{$start}}, $distance."_".$start."_".$stop;
                    }
                    elsif (3 * $start < 3 * length($seq)){
                        next START if ((3*$stop+2 < 3 * length($seq)) && (3*$stop+2 > 3*$start));
                        my $distance = (3*$stop+2) - (3*$start);
                        push @{$hash_start{$start}}, $distance."_".$start."_".$stop;
                    }
                }   
        }
        foreach my $start (@{$frame_startCodon[1]}){
            next unless (exists $hash_start{$start});
            @{$hash_start{$start}} = map {$_ ->[0]} sort {$a->[1] <=> $b->[1]} map {[$_,$_=~/(\d+)\_.+\_.+/]} @{$hash_start{$start}};
#print join"\t",@{$hash_start{$start}},"\n";
            push @{$frame_start_stop[1]}, ${$hash_start{$start}}[0];
        }

        if (!@{$frame_start_stop[1]}){ $frame_noORF = "across0"; }
        else{
            @{$frame_start_stop[1]} = map {$_ ->[0]} sort {$b->[1] <=> $a->[1]} map {[$_,$_=~/(\d+)\_.+\_.+/]} @{$frame_start_stop[1]};
            my @temp = split(/\_/,${$frame_start_stop[1]}[0]);
            my $start = 3 * $temp[1]; my $stop = 3 * $temp[2] + 2;

            if ($start < length($seq)){
                $frame_start = $start;
                if ($stop < 2*length($seq)){
                    $frame_stop = $stop - length($seq);
                    $frame_across = 1;
                }
                elsif ($stop < 3*length($seq)){
                    $frame_stop = $stop - 2*length($seq);
                    $frame_across = 2;
                }
                else{
                    $frame_stop = $stop - 3*length($seq);
                    $frame_across = 1;
                }
            }
            elsif($start < 2*length($seq)){
                $frame_start = $start - length($seq);
                if ($stop < 3*length($seq)){
                    $frame_stop = $stop - 2*length($seq);
                    $frame_across = 1;
                }
                else{
                    $frame_stop = $stop - 3*length($seq);
                    $frame_across = 2;
                }
            }
            elsif($start < 3*length($seq)){
                $frame_start = $start - 2*length($seq);
                $frame_stop = $stop - 3*length($seq);
                $frame_across = 1;
            }

            my $frame_str_a = substr($seq, $frame_start);
            my $frame_str_b = substr($seq, 0, $frame_stop + 1);
            if ($frame_across == 1){ $frame_ORF = $frame_str_a.$frame_str_b;}
            elsif ($frame_across == 2){$frame_ORF = $frame_str_a.$seq.$frame_str_b;}
            
            if(length($frame_ORF) >= 63){
                $frame_noORF = "-";
                return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
            }else{
                $frame_ORF = "";
                $frame_noORF = "less_than20aa";
            }
        }
    }
    
    
    if(!($frame_ORF) && @{$frame_startCodon[2]}){
        my %hash_start;
        START:foreach my $start (@{$frame_startCodon[2]}){
                foreach my $stop (@quaSeq_stopCodon){
                    next if ($stop < $start);
                    if (3 * $start < length($seq)){
                        next START if ((3*$stop+2 < length($seq)) && (3*$stop+2 > 3*$start));
                        my $distance = (3*$stop+2) - (3*$start);
                        push @{$hash_start{$start}}, $distance."_".$start."_".$stop;
                    }
                    elsif (3 * $start < 2 * length($seq)){
                        next START if ((3*$stop+2 < 2 * length($seq)) && (3*$stop+2 > 3*$start));
                        my $distance = (3*$stop+2) - (3*$start);
                        push @{$hash_start{$start}}, $distance."_".$start."_".$stop;
                    }
                    elsif (3 * $start < 3 * length($seq)){
                        next START if ((3*$stop+2 < 3 * length($seq)) && (3*$stop+2 > 3*$start));
                        my $distance = (3*$stop+2) - (3*$start);
                        push @{$hash_start{$start}}, $distance."_".$start."_".$stop;
                    }
                }
        }
        foreach my $start (@{$frame_startCodon[2]}){
            next unless (exists $hash_start{$start});
            @{$hash_start{$start}} = map {$_ ->[0]} sort {$a->[1] <=> $b->[1]} map {[$_,$_=~/(\d+)\_.+\_.+/]} @{$hash_start{$start}};
#print join"\t",@{$hash_start{$start}},"\n";exit;
            push @{$frame_start_stop[2]}, ${$hash_start{$start}}[0];
        }

        if (!@{$frame_start_stop[2]}){ $frame_noORF = "across0";}
        else{
            @{$frame_start_stop[2]} = map {$_ ->[0]} sort {$b->[1] <=> $a->[1]} map {[$_,$_=~/(\d+)\_.+\_.+/]} @{$frame_start_stop[2]};
            my @temp = split(/\_/,${$frame_start_stop[2]}[0]);
            my $start = 3 * $temp[1]; my $stop = 3 * $temp[2] + 2;

            if ($start < length($seq)){
                $frame_start = $start;
                if ($stop < 2*length($seq)){
                    $frame_stop = $stop - length($seq);
                    $frame_across = 1;
                }
                elsif ($stop < 3*length($seq)){
                    $frame_stop = $stop - 2*length($seq);
                    $frame_across = 2;
                }
                else{
                    $frame_stop = $stop - 3*length($seq);
                    $frame_across = 1;
                }
            }
            elsif($start < 2*length($seq)){
                $frame_start = $start - length($seq);
                if ($stop < 3*length($seq)){
                    $frame_stop = $stop - 2*length($seq);
                    $frame_across = 1;
                }
                else{
                    $frame_stop = $stop - 3*length($seq);
                    $frame_across = 2;
                }
            }
            elsif($start < 3*length($seq)){
                $frame_start = $start - 2*length($seq);
                $frame_stop = $stop - 3*length($seq);
                $frame_across = 1;
            }

            my $frame_str_a = substr($seq, $frame_start);
            my $frame_str_b = substr($seq, 0, $frame_stop + 1);
            if ($frame_across == 1){ $frame_ORF = $frame_str_a.$frame_str_b;}
            elsif ($frame_across == 2){$frame_ORF = $frame_str_a.$seq.$frame_str_b;}

            if(length($frame_ORF) >= 63){
                $frame_noORF = "-";
                return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
            }else{
                $frame_ORF = "";
                $frame_noORF = "less_than20aa";
            }
        }
    }


    if(!($frame_ORF) && @{$frame_startCodon[3]}){
        my %hash_start;
        START:foreach my $start (@{$frame_startCodon[3]}){
                foreach my $stop (@quaSeq_stopCodon){
                    next if ($stop < $start);
                    if (3 * $start < length($seq)){
                        next START if ((3*$stop+2 < length($seq)) && (3*$stop+2 > 3*$start));
                        my $distance = (3*$stop+2) - (3*$start);
                        push @{$hash_start{$start}}, $distance."_".$start."_".$stop;
                    }
                    elsif (3 * $start < 2 * length($seq)){
                        next START if ((3*$stop+2 < 2 * length($seq)) && (3*$stop+2 > 3*$start));
                        my $distance = (3*$stop+2) - (3*$start);
                        push @{$hash_start{$start}}, $distance."_".$start."_".$stop;
                    }
                    elsif (3 * $start < 3 * length($seq)){
                        next START if ((3*$stop+2 < 3 * length($seq)) && (3*$stop+2 > 3*$start));
                        my $distance = (3*$stop+2) - (3*$start);
                        push @{$hash_start{$start}}, $distance."_".$start."_".$stop;
                    }
                }
        }
        foreach my $start (@{$frame_startCodon[3]}){
            next unless (exists $hash_start{$start});
            @{$hash_start{$start}} = map {$_ ->[0]} sort {$a->[1] <=> $b->[1]} map {[$_,$_=~/(\d+)\_.+\_.+/]} @{$hash_start{$start}};
            push @{$frame_start_stop[3]}, ${$hash_start{$start}}[0];
        }

        if (!@{$frame_start_stop[3]}){ $frame_noORF = "across0";}
        else{
            @{$frame_start_stop[3]} = map {$_ ->[0]} sort {$b->[1] <=> $a->[1]} map {[$_,$_=~/(\d+)\_.+\_.+/]} @{$frame_start_stop[3]};
            my @temp = split(/\_/,${$frame_start_stop[3]}[0]);
            my $start = 3 * $temp[1]; my $stop = 3 * $temp[2] + 2;

            if ($start < length($seq)){
                $frame_start = $start;
                if ($stop < 2*length($seq)){
                    $frame_stop = $stop - length($seq);
                    $frame_across = 1;
                }
                elsif ($stop < 3*length($seq)){
                    $frame_stop = $stop - 2*length($seq);
                    $frame_across = 2;
                }
                else{
                    $frame_stop = $stop - 3*length($seq); 
                    $frame_across = 1;
                }
            }
            elsif($start < 2*length($seq)){
                $frame_start = $start - length($seq);
                if ($stop < 3*length($seq)){
                    $frame_stop = $stop - 2*length($seq);
                    $frame_across = 1;
                }
                else{
                    $frame_stop = $stop - 3*length($seq);
                    $frame_across = 2;
                }
            }
            elsif($start < 3*length($seq)){
                $frame_start = $start - 2*length($seq);
                $frame_stop = $stop - 3*length($seq);
                $frame_across = 1;
            }

            my $frame_str_a = substr($seq, $frame_start);
            my $frame_str_b = substr($seq, 0, $frame_stop + 1);
            if ($frame_across == 1){ $frame_ORF = $frame_str_a.$frame_str_b;}
            elsif ($frame_across == 2){$frame_ORF = $frame_str_a.$seq.$frame_str_b;}

            if(length($frame_ORF) >= 63){
                $frame_noORF = "-";
                return ($frame_start, $frame_stop, $frame_across, $frame_ORF, $frame_noORF);
            }else{
                $frame_ORF = "";
                $frame_noORF = "less_than20aa";
            }
        }
    }

    if (!($frame_ORF) && (@{$frame_startCodon[0]} || @{$frame_startCodon[1]} || @{$frame_startCodon[2]} || @{$frame_startCodon[3]})){
        return ("-","-","-","-", $frame_noORF);
    }
}

###################################################################################
sub translation{
    my ($ORF) = @_;
    my $peptide = "";
    if ($ORF ne "-"){
        for (my $i=0; $i<=(length($ORF)-3); $i+=3){
            $peptide .= codonMap(substr($ORF,$i,3));
        }
        return $peptide;
    }else {return "-";}
}

##############################################################
sub codonMap{ #input:($codon); output:($aa)
    my %codonMap = (
            TTT => "F", TTC => "F", TTA => "L", TTG => "L",
            TCT => "S", TCC => "S", TCA => "S", TCG => "S",
            TAT => "Y", TAC => "Y", TAA => "*", TAG => "*",
            TGT => "C", TGC => "C", TGA => "*", TGG => "W",
            CTT => "L", CTC => "L", CTA => "L", CTG => "L",
            CCT => "P", CCC => "P", CCA => "P", CCG => "P",
            CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
            CGT => "R", CGC => "R", CGA => "R", CGG => "R",
            ATT => "I", ATC => "I", ATA => "I", ATG => "M",
            ACT => "T", ACC => "T", ACA => "T", ACG => "T",
            AAT => "N", AAC => "N", AAA => "K", AAG => "K",
            AGT => "S", AGC => "S", AGA => "R", AGG => "R",
            GTT => "V", GTC => "V", GTA => "V", GTG => "V",
            GCT => "A", GCC => "A", GCA => "A", GCG => "A",
            GAT => "D", GAC => "D", GAA => "E", GAG => "E",
            GGT => "G", GGC => "G", GGA => "G", GGG => "G",

            TCN => "S", CTN => "L", CCN => "P", CGN => "R",
            ACN => "T", GTN => "V", GCN => "A", GGN => "G",
            );
    my ($codon) = @_;
    $codon = uc $codon; #uppercase
    if(exists $codonMap{$codon}){return $codonMap{$codon};}
    else{
        print BAD "Bad codon \"$codon\"!!\n";
        return "X";
    }
}
### NNN;
##A/C/G/T-NN, N-A/C/G/T-N, NN-A/C/G/T;
#A/C/G/T-A/C/G/T-N, A/C/G/T-N-A/C/G/T, N-A/C/G/T-A/C/G/T;



