#!usr/bin/perl -w
use strict;
print << "EOF";
This program is to split m6A data from REPIC to bed format(non-redundant) by tools.
EOF

die "Usage: perl $0 <m6A_sites_species_human_hg38.txt> <m6A_exomePeak.bed> <m6A_macs2.bed> <m6A_MeTPeak.bed>\n" unless @ARGV==4;

open M6A, "<", "$ARGV[0]" or die ($!);
open EXOP, ">", "$ARGV[1]" or die ($!);
open MACS, ">", "$ARGV[2]" or die ($!);
open METP, ">", "$ARGV[3]" or die ($!);

my %pos;
my $exomePeak = 0;
my $macs2 = 0;
my $MeTPeak = 0;

<M6A>;
while(<M6A>){
    chomp;
    my @line = split(/\t/,$_);
    $line[0] =~ /(?<chrom>chr.+)\:(?<start>\d+)\-(?<end>\d+)\[(?<strand>.+)\]/;
    $pos{"$line[7]_$line[0]"} += 1;
    next if ($pos{"$line[7]_$line[0]"} > 1); #skip duplicate coordinates
    if ($line[7] eq "exomePeak"){
        print EXOP "$+{chrom}\t$+{start}\t$+{end}\t$line[7]\t$line[10]\t$+{strand}\n";
        $exomePeak += 1;
    }
    elsif ($line[7] eq "macs2"){
        print MACS "$+{chrom}\t$+{start}\t$+{end}\t$line[7]\t$line[10]\t$+{strand}\n";
        $macs2 += 1;
    }
    elsif ($line[7] eq "MeTPeak"){
        print METP "$+{chrom}\t$+{start}\t$+{end}\t$line[7]\t$line[10]\t$+{strand}\n";
        $MeTPeak += 1;
    }
}
close M6A;
close EXOP;
close MACS;
close METP;

print "$exomePeak exomePeak\n";
print "$macs2 macs2\n";
print "$MeTPeak MeTPeak\n";
