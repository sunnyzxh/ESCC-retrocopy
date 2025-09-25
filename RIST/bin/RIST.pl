#!/usr/bin/perl -w
use strict;

open IN, "config.txt";
my ($samtools, $fa, $sample_info, $input, $output, $output1, $bam_dir, $refgene);
while (<IN>) {
    chomp;
    my @info = split /=/, $_;
    $samtools = $info[-1] if ( $info[0] eq "samtools" );
    $fa = $info[-1] if ( $info[0] eq "reference" );
    $refgene = $info[-1] if ( $info[0] eq "refgene" );
    $sample_info = $info[-1] if ( $info[0] eq "sample_info" );
    $input = $info[-1] if ( $info[0] eq "gripper_ouput" );
    $output = $info[-1] if ( $info[0] eq "output" );
    $bam_dir = $info[-1] if ( $info[0] eq "bam_dir" );
}
$output1 = "$output\_moredetail.txt";
close IN;

###============Output header=====================###
my $header = "Sample\tID\tParent_Gene(P)\tP_exon_n\tP_chr\tP_exon_start\tP_exon_end\tP_retro_s\tP_retro_e\tRetro_exon_n\tInsert_chr\tInsert_start\tInsert_end\tInsert_depth\tN_read_pairs\tN_junction_reads\tTSDs\tPolyA\tRetro_seq\tConfidence\tPair_Conf\tPair_sample\tPair_ID\tPair_retro_s\tPair_retro_e\tPair_retro_exon_n\tPair_Ins_chr\tPair_Ins_s\tPair_Ins_e\tPair_Ins_depth\tPair_N_read_pairs\tPair_N_junction_reads\tPair_TSDs\tPair_PolyA\tPair_retro_seq";

###============Refgene INFO======================###
open ING, "$refgene";
#785 NM_007276 chr7 + 26241062 26253227 26242618 26251828 6 26241062,26242590,26245987,26248012,26251281,26251701, 26241181,26242642,26246130,26248175,26251376,26253227, 0 CBX3 cmpl cmpl -1,0,0,2,0,2,
my %ref_genes;
while (<ING>) {
    chomp;
    my @info = split;
    next if ( $info[2] =~ /_/ );  ##remove genes in other contigs
    my $line = $_;
    if ( !$ref_genes{$info[12]} ) {
        $ref_genes{$info[12]} = $line;
    }
    else {
        $ref_genes{$info[12]} = "$ref_genes{$info[12]}\n$line"; ##for genes with multiple transcript,need to check if multiple chrs
    }
}
close ING;

###===========Paired Samples=====================###
open INS, "$sample_info";
#R17025773LD01-BDESCCDNA001-T    001     Tumor   baidu2
my %pair_types; $pair_types{"Tumor"} = "Normal"; $pair_types{"Normal"} = "Tumor";
my %patient_type; my %smps; my %id;
while (<INS>) {
    chomp;
    my @info = split;
    $patient_type{"$info[1]\t$info[2]"} = $info[0];
    $smps{$info[0]} = "$info[1]\t$info[2]";
    my $type = "T"; $type = "N" if ( $info[2] eq "Normal" );
    $id{$info[0]} = "$info[1]"."$type";
}
close INS;

###===========Main==============================###
open IN, "$input";
#intronic TMEM132B chr12 125801102 125801147 Germline chr12 104359592 104382656 104359629 104382442 TDG FJ1705100058LD01.sorted.markdup.bqsr
open OUT, ">$output.txt";;
open OUT1, ">$output1";
print OUT1 "$header\n";

while (<IN>) {
    chomp;
    my @info = split;
    my $line = $_;
    ##===Sample info===##
    my $sample = $info[11];
    $info[11] =~ s/.sorted.markdup//g;  ###should be changed according to your file name to make the sampleID the same with that in sample_info
    my ( $patient, $type ) = split /\t/, $smps{$info[11]};
    my $pair0 = $patient_type{"$patient\t$pair_types{$type}"};
    my $pair = "$pair0.sorted.markdup";
    ##===Get bam===##
    my $bam_id = 0; my $bam = ""; my $pairbam = "";
    for my $s ( $sample, $pair ) {
        $bam_id ++;
        my $b =  `ls $bam_dir/$s.bam`; chomp $b;
        print "error $s\n" if ( !-e $b );
        $bam = $b if ( $bam_id == 1 );
        $pairbam = $b if ( $bam_id == 2 );
    }
    ##===Insertion range===##
    my $insert_s = $info[1] - 100;
    my $insert_e = $info[2] + 100;
    my $insert_region  = "$info[0]:$insert_s-$insert_e";
    ##===Parent range===##
    my $parent_s = $info[4] - 100;
    my $parent_e = $info[5] + 100;
    ##===Match chr===##
    my $matchr  = $info[3]; #for inquring insertion sites;
    my $matchr1 = $info[0]; #for inquring parent region;
    if ( $info[0] eq $info[3] ) {
        $matchr = "=";
        $matchr1 = "=";
    }
    ##===Get info===##
    my $count = 0;
    my %exon_n; my %exon_r;
    for my $b ( $bam, $pairbam ) {
        $count ++;
        #==Discordant read pairs==#
        my $dis_rp = 0;
        $dis_rp = `$samtools view -F 1280 $b $insert_region |awk '\$7 == "$matchr" && \$8 > $parent_s && \$8 < $parent_e' |wc -l`; chomp ($dis_rp);
     
        #==Mean Depth==#
        my $mean_depth = 0;
        my $cover = `$samtools coverage $b -r $info[0]:$info[1]-$info[2]`; chomp $cover;
        $mean_depth = ( split /\t/, (split /\n/, $cover)[-1] )[6];

        #==PolyA tail==#
        open INP, "$samtools view -F 1280 $b $info[0]:$info[1]-$info[2] |grep -v 150M |" or die $!;
        my %poly_a; my $poly = "NA";
        while (<INP>) {
            chomp;
            my @info_d = split;
            if ( $info_d[5] =~ /^(\d+)M(\d+)S/ ) {
                my $m = ( $info_d[5] =~ /(\d+)M/ )[0];
                my $sub_s = substr($info_d[9], $m);
                if ( $sub_s =~ /^(T+)[AGC]/ ) {
                    my $ts = ($sub_s =~ /^(T+)[AGC]/)[0];
                    $poly_a{$ts} ++;
                }
            }
        }
        close INP;
        for my $p ( sort keys %poly_a ) {
            if ( length($p) > 3 and $poly_a{$p} > 2 ) { my $l = length($p); $poly = "A($l)"; }
        }
        #==Breakpoints and TSDs==##
        open INF, "$samtools view -F 1280 $b $insert_region |grep SA: |" or die $!;
        my @insert = (); my $ins_s = ""; my $ins_e = ""; my $insertion = "NA"; my $seq = "NA";
        while (<INF>) {
            chomp;
            my @info_d = split;
            if ( $info_d[5] =~ /(\d+)S(\d+)M/ ) {
                push @insert, "$info_d[2],$info_d[3]";
            }
            elsif ( $info_d[5] =~ /(\d+)M(\d+)S/ ) {
                my $m = ( $info_d[5] =~ /(\d+)M(\d+)S/ )[0];
                my $pos = $info_d[3] + $m - 1;
                push @insert, "$info_d[2],$pos";
            }
        }
        close INF;

        open INF, "$samtools view -F 1280 $b $insert_region |" or die $!;
        #ST-E00493:203:H33YLCCXY:8:1215:12449:24251      163     chr15   40854234        60      67S83M  =       40854234        106 seq qua SA:Z:chr4,49138654,+,59M91S,0,3;  PG:Z:MarkDuplicates     AS:i:83 XS:i:20 MD:Z:83 NM:i:0  RG:Z:FP1705100059DN01

#intronic TMEM132B chr12 125801102 125801147 Germline chr12 104359592 104382656 104359629 104382442 TDG FJ1705100058LD01.sorted.markdup.bqsr
        my @parent_s = (); my @parent_e = (); my @parent_e0 = ();
        my $parent_start = ""; my $parent_end = "";
        while (<INF>) {
            chomp;
            my @info_d = split;
            if ( $info_d[5] =~ /(\d+)M(\d+)S/ and $info_d[6] eq $matchr and $info_d[7] > $parent_s and $info_d[7] < $parent_e ) {    
                my $pos_parent_e = $info_d[7];
                push @parent_e0, "$pos_parent_e";
            }
            if ( $info_d[11] =~ /SA:/ ) {
                my @sainfo = split /,/, $info_d[11];
                if ( $sainfo[3] =~ /(\d+)S(\d+)M/ ) {
                    push @parent_s,"$info[3],$sainfo[1]" if ($sainfo[1] > $parent_s and $sainfo[1] < $parent_e and $sainfo[0] eq "SA:Z:$info[3]");
                }
            }
        }
        close INF;
        
        if ( $#parent_e0 == -1 ) {
            $parent_end = $info[5] if ( $count == 1 );
            $parent_end = "NA" if ( $count == 2 );
        }
        else {
            @parent_e0 = sort(@parent_e0);
            open INF1, "$samtools view -F 1280 $b $info[3]:$parent_e0[0]-$parent_e0[-1] |awk '\$7 == \"$matchr1\" && \$8 > $insert_s && \$8 < $insert_e' |" or die $!;
            while (<INF1>) {
                chomp;
                my @info_d = split;
                if ( $info_d[5] =~ /^(\d+)M/ ) {
                    my $m = ( $info_d[5] =~ /(\d+)M/ )[0];
                    my $pos = $info_d[3] + $m - 1;
                    push @parent_e, "$info_d[2],$pos";
                }
            }
            close INF1;
            my $p_e = pos_cluster ( @parent_e, 100 );
            $parent_end = ( split /[,:;]/, $p_e )[1];
        }

        if ( $#parent_s == -1 ) {
            $parent_start = $info[4] if ( $count == 1 );
            $parent_start = "NA" if ( $count == 2 );
        }
        else {
            my $p_s = pos_cluster ( @parent_s, 10 );
            $parent_start = ( split /[,:;]/, $p_s )[1];
        }
        my ( $ins_chr, $ins_s_tmp, $c1_n, $ins_e_tmp, $c2_n ) = ( "NA", 0, 0, 0, 0 );
        if ( $#insert == -1 ) {
            if ( $count == 1 ) {$ins_s = $info[1]; $ins_e = $info[2]; }
            if ( $count == 2 ) {$ins_s = "NA"; $ins_e = "NA"; }
        }
        else {
            $insertion = pos_cluster ( @insert, 5); #chr12,104382489:13;chr12,104382365:2
            ( $ins_chr, $ins_s_tmp, $c1_n, $ins_e_tmp, $c2_n ) = ( split /[,;:]/, $insertion )[0, 1, 2, 4, 5];
            $ins_e_tmp = $ins_s_tmp if ( !$ins_e_tmp ); $c2_n = 0 if ( !$c2_n );
            my @ins_tmp = sort($ins_s_tmp, $ins_e_tmp);
            $ins_s = $ins_tmp[0];
            $ins_e = $ins_tmp[1];
            if ( $ins_e - $ins_s + 1 < 60 ) {
                if ( $fa =~ /demo.fasta/ ) {my $ins_s_demo = $ins_s - 104349629 + 1; my $ins_e_demo = $ins_e - 104349629 + 1; my $ins_chr_demo = "chr12:104349629-125811147"; $seq = `$samtools faidx $fa $ins_chr_demo:$ins_s_demo-$ins_e_demo |tail -1`; chomp $seq; }
                else { $seq = `$samtools faidx $fa $ins_chr:$ins_s-$ins_e |tail -1`; chomp $seq; }
            }
        }
       
        ##==Exon-exon junction reads and Retro sequences==##
        #785 NM_007276 chr7 + 26241062 26253227 26242618 26251828 6 26241062,26242590,26245987,26248012,26251281,26251701, 26241181,26242642,26246130,26248175,26251376,26253227, 0 CBX3 cmpl cmpl -1,0,0,2,0,2,
        #intronic TMEM132B chr12 125801102 125801147 Germline chr12 104359592 104382656 104359629 104382442 TDG FJ1705100058LD01.sorted.markdup.bqsr
        my $n_junction = 0; my $retro_seq = ""; my $transcript = ""; my $trans_exon_n = 1;
        if ( !$ref_genes{$info[7]} ) {
            $n_junction = "NA";
            my $seq0;
            if ( $fa !~ "demo.fasta" ) {$seq0 = `$samtools faidx $fa $info[3]:$info[4]-$info[5]`; chomp $seq0;}
            if ( $fa =~ /demo.fasta/ ) {my $demo_c = "chr12:104349629-125811147"; my $demo_s = $info[4]-104349629 + 1; my $demo_e = $info[5]-104349629 + 1; $seq0 = `$samtools faidx $fa $demo_c:$demo_s-$demo_e`; chomp $seq0;}
            my @seqs = split /\n/, $seq0;
            for my $s ( 1..$#seqs ) {
                $retro_seq = "$retro_seq"."$seqs[$s]";
            }
            $exon_n{$info[7]} = 1;
            $exon_r{$info[7]} = "$info[3]\t$info[4]\t$info[5]";
        }
        elsif ( $parent_start ne "NA" and $parent_end ne "NA" ) {
            my @trans = split /\n/, $ref_genes{$info[7]};
            my $n = 0;
            for my $t ( @trans ) {
                my @info_t = split /\t/, $t;
                next if ( $info_t[2] ne $info[3] );
                open INJ, "$samtools view -F 1280 $b $info_t[2]:$info_t[4]-$info_t[5] |awk '\$6 != \"150M\"' |" or die $!; 
                while (<INJ>) {
                    chomp;
                    my @info_d = split;
                    my $pos = ""; my $flag = 0;
                    if ( $info_d[5] =~ /^(\d+)S(\d+)M/ ) {
                        $pos = $info_d[3];
                    }
                    elsif ( $info_d[5] =~ /^(\d+)M(\d+)S/ ) {
                        my $m = ( $info_d[5] =~ /(\d+)M(\d+)S/ )[0];
                        $pos = $info_d[3] + $m - 1;
                    } 
                    else { next; }
                    for my $p ( (split /,/, $info_t[9]), (split /,/, $info_t[10]) ) {
                        $flag = 1 if ( abs($pos-$p) < 10 );
                    }
                    $n ++ if ( $flag == 1 );
                }
                close INJ;
                if ( $n_junction <= $n ) {
                    $n_junction = $n;
                    $transcript = $t;
                }
            }
            #=get retro sequences=#
            my @trans_d = split /\t/, $transcript;
            $exon_n{$trans_d[12]} = $trans_d[8];
            $exon_r{$trans_d[12]} = "$trans_d[2]\t$trans_d[4]\t$trans_d[5]";
            my @exon_s = split /,/, $trans_d[9];
            my @exon_e = split /,/, $trans_d[10];
            my ( $start_exon, $end_exon ) = (0, 0);
            for my $i ( 0..$#exon_s ) {
                $start_exon = $i if ( abs($exon_s[$i] - $parent_start) < 10 );
                $end_exon = $i if ( $parent_end >= $exon_s[$i] and $parent_end <= $exon_e[$i] );
            }
            if ( $parent_start < $exon_s[0] ) { $start_exon = 0; $parent_start = $exon_s[0]; }
            if ( $parent_end > $exon_e[-1] ) { $end_exon = $#exon_s; $parent_end = $exon_e[-1]; }
            $trans_exon_n = $end_exon - $start_exon + 1;
            for my $i ( $start_exon..$end_exon ) {
                my $seq0;
                if ( $fa !~ "demo.fasta" ) {$seq0 = `$samtools faidx $fa $trans_d[2]:$exon_s[$i]-$exon_e[$i]`; chomp $seq0;}
                if ( $fa =~ /demo.fasta/ ) {my $demo_c = "chr12:104349629-125811147"; my $demo_s = $exon_s[$i]-104349629 + 1; my $demo_e = $exon_e[$i]-104349629 + 1; $seq0 = `$samtools faidx $fa $demo_c:$demo_s-$demo_e`; chomp $seq0;}
                my @seqs = split /\n/, $seq0;
                for my $s ( 1..$#seqs ) {
                    $retro_seq = "$retro_seq"."$seqs[$s]";
                }
            }
        }
        ##==decide whether sample have this retroduplication==##
        my $pairinfo = "pairNO";
        $pairinfo = "pairYES" if ( $dis_rp > 0 and ($c1_n > 1 or $c2_n > 1) and ($dis_rp < 1000 and $c1_n < 1000 and $c2_n < 1000) );
        $pairinfo = "not_sure" if ( $mean_depth < 5 );
        my $tinfo = "YES";
        $tinfo = "NO" if ( $dis_rp > 1000 or $c1_n > 1000 or $c2_n > 1000 );
        ##==OUTPUT==##
        $n_junction = "NA" if ( $trans_exon_n == 1 );
        if ( $count == 1 ) {
            print OUT "$line\t$tinfo\t$dis_rp\t$insertion\t$ins_s\t$ins_e\t$mean_depth\t$seq";
            print OUT1 "$info[11]\t$id{$info[11]}\t$info[7]\t$exon_n{$info[7]}\t$exon_r{$info[7]}\t$parent_start\t$parent_end\t$trans_exon_n\t$info[0]\t$ins_s\t$ins_e\t$mean_depth\t$dis_rp\t$n_junction\t$seq\t$poly\t$retro_seq\t$tinfo";
        }
        if ( $count == 2 ) {
            print OUT "\t$pairinfo\t$dis_rp\t$insertion\t$ins_s\t$ins_e\t$mean_depth\t$seq\t$pair\n";
            print OUT1 "\t$pairinfo\t$pair0\t$id{$pair0}\t$parent_start\t$parent_end\t$trans_exon_n\t$info[0]\t$ins_s\t$ins_e\t$mean_depth\t$dis_rp\t$n_junction\t$seq\t$poly\t$retro_seq\n";
        }
    }
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
