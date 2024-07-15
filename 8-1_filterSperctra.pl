#!usr/bin/perl -w
use strict;
print << "EOF";
This program is to filter spectra result from pFind.
1.q-value <= 0.01
2.missed cleavage sites ≤3
3.peptides length ≥8 with a new sequence of at least 2 aa at either side of the back splice junction
4.allowing only common modifications: 
    Carbamidomethyl[C]: cysteine carbamidomethylation
    Oxidation[M]: oxidation of methionine
    Acetyl[ProteinN-term]: protein N-terminal acetylation
    Deamidated[N]: deamidation of asparagine (Asn, N) residues
    Phospho[S], Phospho[T], Phospho[Y]: phosphorylation of serine, threonine, and tyrosine residues.

EOF
die "Usage: perl $0 <circRNA.seq> <pFind.spectra> <out.spectra> \n" unless @ARGV==3;

open CIRC, "<", "$ARGV[0]" or die ($!);
open SPEC, "<", "$ARGV[1]" or die ($!);
open OUT, ">", "$ARGV[2]" or die ($!);

my @common_mods = qw(Carbamidomethyl[C] Oxidation[M] Acetyl[ProteinN-term] Deamidated[N] Phospho[S] Phospho[T] Phospho[Y]);
my %hash_common_mods = map{$_=>1} @common_mods;

###############get sequence of each circRNA
my %seq;
while(<CIRC>){
    chomp;
    my @line = split(/\t/,$_);
    next if ($line[13] eq "-");
    $seq{"$line[3]_frame0"}[0]  = length($line[13]); #key = circID_frame, value[0] = length(seq)
    $seq{"$line[3]_frame1"}[0]  = length($line[13]);
    $seq{"$line[3]_frame2"}[0]  = length($line[13]);
    $seq{"$line[3]_frame0"}[1] = $line[16]; #key = circID_frame, value[1] = across
    $seq{"$line[3]_frame1"}[1] = $line[21];
    $seq{"$line[3]_frame2"}[1] = $line[26];
    $seq{"$line[3]_frame0"}[2] = $line[18]; #key = circID_frame, value[2] = peptide
    $seq{"$line[3]_frame1"}[2] = $line[23];
    $seq{"$line[3]_frame2"}[2] = $line[28];

}
close CIRC;

#####################filter
my $title = <SPEC>;
print OUT "filt_proteins\t$title";

my %circ;
SPECTRA:while(<SPEC>){
    chomp;
    my @line = split(/\t/,$_);
    ###1,2,3-1
    next if ($line[4] > 0.01 || length($line[5]) < 8 || $line[16] > 3);
    ###4
    my @Modifications;
    if ($line[10]){
        my @mods = split(/\;/, $line[10]);
        foreach my $pos_mod (@mods){
            my @temp = split(/\,/, $pos_mod);
            push @Modifications, $temp[1];
        }
        my @diffMod = grep {!$hash_common_mods{$_}} @Modifications;
        next SPECTRA if (@diffMod);
    }
    ###3-2
    my @proteins = split(/\//, $line[12]);
    my @filter_proteins;
    foreach my $pro (@proteins){
        $pro =~ /(?<hg38>.+)\_(?<circID>hsa.+)\_(?<frame>.+)\:(?<start>\d+)\-(?<end>\d+)/;
        my @peps;
        if($seq{"$+{circID}_$+{frame}"}[1] >= 1){
            my $bsj = int(($seq{"$+{circID}_$+{frame}"}[0] -1 - $+{start}) / 3);
            if ($bsj == 0){
                my $pep = substr($seq{"$+{circID}_$+{frame}"}[2], 0, 2);
                push @peps, $pep;
            }
            elsif($bsj == length($seq{"$+{circID}_$+{frame}"}[2])-1){
                my $pep = substr($seq{"$+{circID}_$+{frame}"}[2], length($seq{"$+{circID}_$+{frame}"}[2])-2, 2);
                push @peps, $pep;
            }else{
                my $pep = substr($seq{"$+{circID}_$+{frame}"}[2], $bsj-1, 3);
                push @peps, $pep;
            }
        }
        if($seq{"$+{circID}_$+{frame}"}[1] >= 2){
            my $bsj = int((2*$seq{"$+{circID}_$+{frame}"}[0] -1 - $+{start}) / 3);
            if ($bsj == 0){
                my $pep = substr($seq{"$+{circID}_$+{frame}"}[2], 0, 2);
                push @peps, $pep;
            }
            elsif($bsj == length($seq{"$+{circID}_$+{frame}"}[2])-1){
                my $pep = substr($seq{"$+{circID}_$+{frame}"}[2], length($seq{"$+{circID}_$+{frame}"}[2])-2, 2);
                push @peps, $pep;
            }
            else{
                my $pep = substr($seq{"$+{circID}_$+{frame}"}[2], $bsj-1, 3);
                push @peps, $pep;
            }
        }
        if($seq{"$+{circID}_$+{frame}"}[1] == 3){
            my $bsj = int((3*$seq{"$+{circID}_$+{frame}"}[0] -1 - $+{start}) / 3);
            if ($bsj == 0){
                my $pep = substr($seq{"$+{circID}_$+{frame}"}[2], 0, 2);
                push @peps, $pep;
            }
            elsif($bsj == length($seq{"$+{circID}_$+{frame}"}[2])-1){
                my $pep = substr($seq{"$+{circID}_$+{frame}"}[2], length($seq{"$+{circID}_$+{frame}"}[2])-2, 2);
                push @peps, $pep;
            }
            else{
                my $pep = substr($seq{"$+{circID}_$+{frame}"}[2], $bsj-1, 3);
                push @peps, $pep;
            }
        }
        my $flag = 0;
        foreach my $pep (@peps){
            $flag = 1 if ($line[5] =~ /$pep/);
        }
        push @filter_proteins, $pro if ($flag);
        $circ{$+{circID}} += 1 if($flag);
    }
    next unless(@filter_proteins);
# $circ{$+{circID}} += 1;
    my $unshift = join"\/",@filter_proteins;
    print OUT join"\t",($unshift, @line),"\n";
}
close OUT;
close SPEC;

my $circ = keys %circ;
print "$circ circRNA(s) have MS support\n";

