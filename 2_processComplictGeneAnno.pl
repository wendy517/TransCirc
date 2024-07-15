#usr/bin/perl -w
use strict;
print << "EOF";
This program is to cope with comflict gene annatation for circRNA in circAtlas(without consideration for strand). 
EOF
die "Usage: perl $0  <in> <out> \n" unless @ARGV==2;

open IN, "<", "$ARGV[0]" or die ($!);
open OUT, ">", "$ARGV[1]" or die ($!);

my %circRNAs;
#my $full;
while(<IN>){
    chomp;
    my @line = split(/\t/,$_);
    next if ($line[14] eq "partial"); #425070
    $circRNAs{$line[3]}[0] += 1; #count
    push @{$circRNAs{$line[3]}[1]}, $line[18]; #gene_name
    push @{$circRNAs{$line[3]}[2]}, $line[19]; #gene_id
    $circRNAs{$line[3]}[5] = $line[0]; #chr
    $circRNAs{$line[3]}[6] = $line[1]; #start
    $circRNAs{$line[3]}[7] = $line[2]; #end
    $circRNAs{$line[3]}[8] = $line[5]; #strand
}
close IN;
#print $full,"\n";exit;

open IN, "<", "$ARGV[0]" or die ($!);
while(<IN>){
    chomp;
    my @line = split(/\t/,$_);
    next unless (exists $circRNAs{$line[3]});
    $line[3] =~ /hsa\-(?<gene_symbol>.+)\_\d+/;
    if ($circRNAs{$line[3]}[0] == 1){ #regardless of compliction
        $circRNAs{$line[3]}[3] = ${$circRNAs{$line[3]}[1]}[0]; #gene_name
        $circRNAs{$line[3]}[4] = ${$circRNAs{$line[3]}[2]}[0]; #gene_id
    }
    elsif ($circRNAs{$line[3]}[0] > 1){ #redundant gene annotation
        if ( grep {$_ =~ /$+{gene_symbol}/} @{$circRNAs{$line[3]}[1]}){
            for (my $i = 0; $i < @{$circRNAs{$line[3]}[1]}; $i++) {
                if (${$circRNAs{$line[3]}[1]}[$i] =~ /$+{gene_symbol}/){
                    $circRNAs{$line[3]}[3] = ${$circRNAs{$line[3]}[1]}[$i];
                    $circRNAs{$line[3]}[4] = ${$circRNAs{$line[3]}[2]}[$i];
                }
            }
        }else{
            $circRNAs{$line[3]}[3] = join"\,",@{$circRNAs{$line[3]}[1]};
            $circRNAs{$line[3]}[4] = join"\,",@{$circRNAs{$line[3]}[2]};
        } 
    }
}
close IN;
#print join"\,",@{$circRNAs{"hsa-A1BG-AS1_0001"}[1]},"\n"; print $circRNAs{"hsa-A1BG-AS1_0001"}[3],"\t",$circRNAs{"hsa-A1BG-AS1_0001"}[4],"\n";exit;

my @circRNA = keys %circRNAs; #408882
foreach my $circRNA(@circRNA){
    print OUT "$circRNAs{$circRNA}[5]\t$circRNAs{$circRNA}[6]\t$circRNAs{$circRNA}[7]\t$circRNA\t.\t$circRNAs{$circRNA}[8]\t$circRNAs{$circRNA}[3]\t$circRNAs{$circRNA}[4]\n";
}
close OUT;
