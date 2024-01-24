#!/usr/bin/perl -w
use strict;

open IN, "config.txt";
my ($samtools, $fa, $sample_info, $gripper_output, $output, $bam_dir);
while (<IN>) {
    chomp;
    my @info = split /=/, $_;
    $samtools = $info[-1] if ( $info[0] eq "samtools" );
    $fa = $info[-1] if ( $info[0] eq "reference" );
    $sample_info = $info[-1] if ( $info[0] eq "sample_info" );
    $gripper_output = $info[-1] if ( $info[0] eq "gripper_ouput" );
    $output = $info[-1] if ( $info[0] eq "output" );
    $bam_dir = $info[-1] if ( $info[0] eq "bam_dir" );
}
close IN;

open INS, "$sample_info";
my %types; $types{"Tumor"} = "Normal"; $types{"Normal"} = "Tumor";
my %patient_type; my %smps;
while (<INS>) {
    chomp;
    my @info = split;
    $patient_type{"$info[1]\t$info[2]"} = $info[0];
    $smps{$info[0]} = "$info[1]\t$info[2]";
}
close INS;

open IN, "$gripper_output";
open OUT, ">$output";
while (<IN>) {
    chomp;
    my $line = $_;
    print OUT "$line";
    my @info = split;
    $info[11] =~ s/.sorted.markdup//g;  ###should be changed according to your file name to make the sampleID the same with that in sample_info
    my $sample = $info[11];
    my $insert_s = $info[1] - 100;
    my $insert_e = $info[2] + 100;
    my $insert_region  = "$info[0]:$insert_s-$insert_e";   # e.g. C793C.ARF3-NM_001659
    my @e_dis = ( abs($info[1]-$info[4]), abs($info[1]-$info[5]), abs($info[2]-$info[4]), abs($info[2]-$info[5]) );
    @e_dis = sort @e_dis;
    my $left_most = $e_dis[0] - 1000;
    my $right_most = $e_dis[-1] + 1000;    
    
    my $bam = `ls $bam_dir/$sample*.bam`;
    chomp $bam;
    #if ( !-e $bam ) { $bam = `ls /lustre/home/ICR2/database/bd_cram/work/*/$sample*cram`; chomp $bam; print "error $sample\n" if ( !-e $bam ); }
    my ( $patient, $type ) = split /\t/, $smps{$info[11]};
    my $pair = $patient_type{"$patient\t$types{$type}"};
    print "$pair\n";
    my $pairbam = `ls $bam_dir/$pair*.bam`;    
    chomp $pairbam;
    #if ( !-e $pairbam ) { $pairbam = `ls /lustre/home/ICR2/database/bd_cram/work/*/$pair*cram`; chomp $pairbam; print "error $pairbam\n" if ( !-e $pairbam ); }

    my $count = 0;
    my $output = "NA\tNA\tNA\tNA\tNA"; my $pairinfo = "pairNO"; my $tinfo = "YES";
    for my $b ( $bam, $pairbam ) {
        $count ++;
        my $dis_rp;
        if ( $info[0] ne $info[3] ) {
            $dis_rp = `$samtools view -F 1280 $b $insert_region |awk '\$7 == "$info[3]"' |wc -l`; chomp $dis_rp;
        }
        elsif ( $info[0] eq $info[3] ) {
            $dis_rp = `$samtools view -F 1280 $b $insert_region |awk 'function abs(v) {return v < 0 ? -v : v} abs(\$9)> $left_most && abs(\$9) < $right_most' |wc -l`; chomp $dis_rp;
        }
        open INF, "$samtools view -F 1280 $b $insert_region |grep SA: |" or die $!;
        my @insert = ();
        while (<INF>) {
            chomp;
            my @info1 = split;
            #next if ( ($info[2] ne $info[6]) and ($info1[6] ne $info[6]) and ($info1[11] !~ /$info[6]/) );
            #next if ( ($info[2] eq $info[6]) and (abs($info[8]) < $left_most or abs($info[8]) > $right_most) );
            if ( $info1[5] =~ /(\d+)S(\d+)M/ ) {
                push @insert, "$info1[2],$info1[3]";
            }
            elsif ( $info1[5] =~ /(\d+)M(\d+)S/ ) {
                my $m = ( $info1[5] =~ /(\d+)M(\d+)S/ )[0];
                my $pos = $info1[3] + $m - 1;
                push @insert, "$info1[2],$pos";
            }
        }
        close INF;
        if ( $#insert == -1 ) { 
            print OUT "\t$tinfo\t$dis_rp\tNA\tNA\tNA\tNA" if ( $count == 1 );
            print OUT "\t$pairinfo\t$dis_rp\tNA\tNA\tNA\tNA" if ( $count == 2 );
            next;
        }

        my $insertion = pos_cluster ( @insert, 5);
        my ( $ins_chr, $ins_s_tmp, $c1_n, $ins_e_tmp, $c2_n ) = ( split /[,;:]/, $insertion )[0, 1, 2, 4, 5];
        $c2_n = 0 if ( !$c2_n );
        $ins_e_tmp = $ins_s_tmp if ( !$ins_e_tmp );
        my $ins_s = $ins_s_tmp;
        my $ins_e = $ins_e_tmp;
        if ( $ins_s_tmp > $ins_e_tmp ) { $ins_s = $ins_e_tmp; $ins_e = $ins_s_tmp; }
        my $seq = `$samtools faidx $fa $ins_chr:$ins_s-$ins_e |tail -1`;
        chomp $seq;
        $seq = "NA" if ( $ins_e - $ins_s + 1 >= 60 );
        $output = "$dis_rp\t$insertion\t$ins_s\t$ins_e\t$seq";   

        if ( $dis_rp > 0 and ($c1_n > 1 or $c2_n > 1) and ($dis_rp < 1000 and $c1_n < 1000 and $c2_n < 1000) ) {
            $pairinfo = "pairYES" if ( $count == 2 );
        }
        if ( $dis_rp > 1000 or $c1_n > 1000 or $c2_n > 1000 ) {
            $tinfo = "NO" if ( $count == 1 );
        }
        print OUT "\t$pairinfo\t$output" if ( $count == 2 );
        print OUT "\t$tinfo\t$output" if ( $count == 1 );
    }
    print OUT "\t$pair.sorted.markdup\n";
}

sub pos_cluster {
    my ( @info ) = @_;
    my %hash;
    $hash{$info[0]} = 1;
    for ( my $i = 1; $i <= $#info-1; $i ++ ) {
        my ( $chr, $pos ) = ( split /,/, $info[$i] );
        my $flag = 0;
        for my $k ( sort keys %hash ) {
            my ( $chr_o, $pos_o ) = ( split /,/, $k )[0, 1];
            if ( $chr eq $chr_o and ( abs ($pos - $pos_o) <= $info[-1] )) {
                $hash{$k} ++;
                $flag = 1;
            }
        }
        $hash{$info[$i]} = 1 if ( $flag == 0 );
    }
    my @merge = ();
    my $count = 0;
    for my $k ( sort { $hash{$b} <=> $hash{$a} } keys %hash ) {
        $count ++;
        if ( $count <= 3 ) {
            push @merge, "$k:$hash{$k}";
        }
    }
    my $merges = join ";", @merge;
    $merges = 0 if ( !$merges );
    return ( $merges );
}
