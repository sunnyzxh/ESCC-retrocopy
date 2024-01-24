#!/usr/bin/perl -w
use strict;
use List::MoreUtils qw(any uniq);

my %patients;
open IN, "sample_info.txt";
while (<IN>) {
    chomp;
    my @info = split;
    $patients{$info[0]} = $info[1];
}
close IN;

open IN, "S3_Table.txt";
open OUT, ">S2_Table.txt";

my %samples; my %insertion;
while (<IN>) {
    chomp;
    my @info = split;
    $info[12] =~ s/.sorted.markdup.bqsr//g;
    if ( !$samples{$info[11]}{$info[5]} ) { $samples{$info[11]}{$info[5]} = $info[12]; }
    else { $samples{$info[11]}{$info[5]} = "$samples{$info[11]}{$info[5]},$info[12]"; }

    if ( !$insertion{$info[11]}{$info[5]} ) { $insertion{$info[11]}{$info[5]} = "$info[2]-$info[19]-$info[20]"; }
    else { $insertion{$info[11]}{$info[5]} = "$insertion{$info[11]}{$info[5]},$info[2]-$info[19]-$info[20]"; }
}
close IN;
print OUT "Gene\tTotal_retrocopies\tSamples_count\tSample_frq\tPatients_count\tPatients_frq\tGermline_count\tGermline_sample\tSomatic_count\tSomatic_sample\tLoss_count\tLoss_sample\tGermline_count_distribute\tSomatic_count_distribute\tLoss_count_distribute\tGermline_insertion_cluster\tSomatic_insertion_cluster\tLoss_insertion_cluster\n";

for my $g ( sort keys %samples ) {
    my $total_retro = 0;
    for my $t ( "Germline", "Somatic", "Loss" ) {
        $samples{$g}{$t} = "nosample" if ( !$samples{$g}{$t} );
        $insertion{$g}{$t} = "NA" if ( !$insertion{$g}{$t} );
        my $size = scalar ( grep(!/nosample/, (split /,/, $samples{$g}{$t})) );
        $total_retro += $size;
    }
    my @germline = grep(!/nosample/, (split /,/, $samples{$g}{"Germline"})); my $n_germline = scalar@germline; my $n_germline_uniq = scalar(uniq@germline);
    my @somatic = grep(!/nosample/, (split /,/, $samples{$g}{"Somatic"})); my $n_somatic = scalar@somatic; my $n_somatic_uniq = scalar(uniq@somatic);
    my @loss = grep(!/nosample/, (split /,/, $samples{$g}{"Loss"})); my $n_loss = scalar@loss; my $n_loss_uniq = scalar(uniq@loss);

    my $germline_cluster = sample_cluster($samples{$g}{"Germline"});
    my $somatic_cluster = sample_cluster($samples{$g}{"Somatic"});
    my $loss_cluster = sample_cluster($samples{$g}{"Loss"});

    my $germline_pos_cluster = pos_cluster($insertion{$g}{"Germline"});
    my $somatic_pos_cluster = pos_cluster($insertion{$g}{"Somatic"});
    my $loss_pos_cluster = pos_cluster($insertion{$g}{"Loss"});

    my @total_samples = (@germline, @somatic, @loss);
    my @total_patients = ();
    for my $s ( @total_samples ) {
        push @total_patients, $patients{$s};
    }
    my $total_sample = scalar(uniq@total_samples);  
    my $total_sample_frq = $total_sample / 1326; 
    $total_sample_frq = sprintf("%.4f",$total_sample_frq);
    my $total_patient = scalar(uniq@total_patients); 
    my $total_patient_frq = $total_patient / 663;
    $total_patient_frq = sprintf("%.4f", $total_patient_frq);
    
    print OUT "$g\t$total_retro\t$total_sample\t$total_sample_frq\t$total_patient\t$total_patient_frq\t$n_germline\t$n_germline_uniq\t$n_somatic\t$n_somatic_uniq\t$n_loss\t$n_loss_uniq\t$germline_cluster\t$somatic_cluster\t$loss_cluster\t$germline_pos_cluster\t$somatic_pos_cluster\t$loss_pos_cluster\n";
}

sub pos_cluster {
    my (@info) = @_;
    my %hash;
    my @info1 = split /,/, $info[0];
    if ( $info1[0] eq "NA" ) {
        my $r = "NA";
        return($r);
    }
    @info1 = sort(@info1);
    $hash{ $info1[0] } = 1;
    for ( my $i = 1 ; $i <= $#info1 ; $i++ ) {
        my ( $chr, $pos1, $pos2 ) = ( split /-/, $info1[$i] );
        my $flag = 0;
        for my $k ( sort keys %hash ) {
            my ( $chr_o, $pos_o1, $pos_o2 ) = ( split /-/, $k )[ 0, 1, 2 ];
            if ( $chr eq $chr_o and ( abs( $pos1 - $pos_o1 ) <= 200 ) ) {
                $hash{$k}++;
                $flag = 1;
            }
        }
        $hash{ $info1[$i] } = 1 if ( $flag == 0 );
    }
    my @merge = ();
    for my $k ( sort { $hash{$b} <=> $hash{$a} } keys %hash ) {
        push @merge, "$k,$hash{$k}";
    }

    my $m = join ";", @merge;
    return($m);
}

sub sample_cluster {
    my (@info) = @_;
    my %hash;
    my @info1 = split /,/, $info[0];
    if ( $info1[0] eq "nosample" ) {
        my $r = "NA";
        return($r);
    }
    for my $s ( @info1 ) {
        $hash{$s} ++;
    }
    my %count;
    for my $s ( sort keys %hash ) {
        $count{$hash{$s}} ++;
    }
    my @merge;
    for my $s ( sort { $a <=> $b } keys %count ) {
        push @merge, "$s:$count{$s}";
    }
    my $m = join ";", @merge;
    return($m);
}
