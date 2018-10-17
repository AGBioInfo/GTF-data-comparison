#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use DBI;
use File::Basename;
use Data::Dumper;

my %options;
my $results = GetOptions( \%options,
                          'input_gtf1|g1=s',
                           'input_gtf2|g2=s');
			  #'output|o=s');

####Author Ankur Ganveer##########
################# GLOBALS ##################################
my $input_gtf1;
my $input_gtf2;
my $converted_gtf1 ="intern1.gtf";
my $converted_gtf2 ="intern2.gtf";
my $unique_gtf1 ="intern1_uniq.gtf";
my $unique_gtf2 ="intern2_uniq.gtf";
my $intersect="intersect.gtf";
my $nover_intr="noverlap.gtf";
my $over_intr="overlap.gtf";
my $over_par="overlap_par.gtf";
my $over_par_uniq="overlap_par_uniq.gtf";
my $over_par_diff="overlap_par_diff.gtf";
my $over_diff1 ="overlap_diff1.gtf";
my $over_diff2 ="overlap_diff2.gtf";
my $over_gene="overlap_gene.gtf";
my $over_tran="overlap_trans.gtf";
my $over_exon="overlap_exon.gtf";
my $outputfile1 = "noverlap_count.txt";
my $outputfile2 = "overlap_count.txt";
my $file1_chr = "chr_file1.txt";
my $file2_chr = "chr_file2.txt";
my $file1_snd = "snd_file1.txt";
my $file2_snd = "snd_file2.txt";
my $gtffile1 = "gtf1_count.txt";
my $gtffile2 = "gtf2_count.txt";
my $gtf1_chr = "chr_gtf1.txt";
my $gtf2_chr = "chr_gtf2.txt";
my $gtf1_snd = "snd_gtf1.txt";
my $gtf2_snd = "snd_gtf2.txt";
my $chrfile1 = "chr_intern1.txt";
my $chrfile2 = "chr_intern2.txt";
############################################################

################# VERIFYING COMMAND OPTIONS ################
&check_options( \%options );
############################################################
######Opening the files and getting the data################
############################################################

### Check and convert input files in a processing format ###
my $conv_gtf1 = convert($input_gtf1,$converted_gtf1);
my $conv_gtf2 = convert($input_gtf2,$converted_gtf2);

### Remove duplicates from <fileA> <fileB> ###
my $uniq_gtf1 = `awk '!x[\$0]++' $conv_gtf1 >$unique_gtf1`;
my $uniq_gtf2 = `awk '!x[\$0]++' $conv_gtf2 >$unique_gtf2`;

### Find all the unique chromosomes from <fileA> and <fileB> ###
my $chr1 =`awk < $unique_gtf1 '{print \$1}' | sort | uniq -c > $chrfile1`;
my $chr2 =`awk < $unique_gtf2 '{print \$1}' | sort | uniq -c > $chrfile2`;

### Find overlapping features of <fileA> in <fileB> using intersectBed tool of BedTools ###
my $lines = `intersectBed -a $unique_gtf1 -b $unique_gtf2 -wao -f 0.50 -r -s >$intersect`;

### Differentiating between overlapping and non-overlapping features ###
my $nover = `awk '{if(\$10==".") print\$0}' $intersect > $nover_intr`;
my $over = `awk '{if(\$10!=".") print\$0}' $intersect > $over_intr`;

### Reducing overlapped features to only same feature type ###
my $overpar = `awk '{if(\$3==\$12) print\$0}' $over_intr > $over_par`;

### Removing the repitition of an entry based on <fileA> entry ###
my $overparuniq = `awk '!x[\$1,\$3,\$4,\$5,\$7]++' $over_par > $over_par_uniq`;

### Finding the differences of corresponding start and corresponding end positions of overlapped feature ###
my $swadif = swap_diff($over_par_uniq,$over_par_diff);

### Creating separate files for gene, transcript and exon features for histogram plots
### Differences if restricted to minimum of -1000 and maximum of 1000
my $ovrdiff1 = `awk '{if((\$14<=1000)&&(\$14>=-1000)) print\$0}' $over_par_diff > $over_diff1`;
my $ovrdiff2 = `awk '{if((\$15<=1000)&&(\$15>=-1000)) print\$0}' $over_diff1 > $over_diff2`;
my $ovrgene = `awk '{if(\$3=="gene") print\$0}' $over_diff2 > $over_gene`;
my $ovrtran = `awk '{if(\$3=="transcript") print\$0}' $over_diff2 > $over_tran`;
my $ovrexon = `awk '{if(\$3=="exon") print\$0}' $over_diff2 > $over_exon`;

### Finding the count of different features in overlap and non-overlap
my ($filecount1,$opr1)=count($nover_intr);
my ($filecount2,$opr2)=count($over_par_diff);

### Printing the count to tab-deliminated text file
feafile($outputfile1,$filecount1,$opr1);
feafile($outputfile2,$filecount2,$opr2);

### Breaking the count of different features by chromosome type in overlap and non-overlap
my ($filechrcount1,$opr_chr1)=chrcount($nover_intr);
my ($filechrcount2,$opr_chr2)=chrcount($over_par_diff);

### Breaking the count of different features by strand type in overlap and non-overlap
my ($filesndcount1,$opr_snd1)=sndcount($nover_intr);
my ($filesndcount2,$opr_snd2)=sndcount($over_par_diff);

### Finding the count of different features in <fileA> and <fileB> ###
my ($gtfcount1,$gtf1)=count($unique_gtf1);
my ($gtfcount2,$gtf2)=count($unique_gtf2);

### Printing the count to tab-deliminated text file ###
feafile($gtffile1,$gtfcount1,$gtf1);
feafile($gtffile2,$gtfcount2,$gtf2);

### Breaking the count of different features by chromosome type in <fileA> and <fileB> ###
my ($gtfchrcount1,$gtf_chr1)=chrcount($unique_gtf1);
my ($gtfchrcount2,$gtf_chr2)=chrcount($unique_gtf2);

### Breaking the count of different features by strand type in <fileA> and <fileB> ###
my ($gtfsndcount1,$gtf_snd1)=sndcount($unique_gtf1);
my ($gtfsndcount2,$gtf_snd2)=sndcount($unique_gtf2);

### Creating list of hashes of breakdown counts ###
my @hash=($opr_chr1,$opr_chr2,$opr_snd1,$opr_snd2,$gtf_chr1,$gtf_chr2,$gtf_snd1,$gtf_snd2);

### Creating list of files ###
my @fh=($file1_chr,$file2_chr,$file1_snd,$file2_snd,$gtf1_chr,$gtf2_chr,$gtf1_snd,$gtf2_snd);

### Printing the breakdown counts in separate files ###
for (my $i=0; $i<=7; $i++){
    my $file=$fh[$i];
    open our ($f), ">$file";
    tabout($hash[$i]);
    close $f;
    }

############## R-Script is been called here ###########################
my $rs = `Rscript compare.R $input_gtf1 $input_gtf2`;
#######################################################################

####################DELETES all the files created during the processing###########################
system("rm intern1.gtf intern2.gtf intern1_uniq.gtf intern2_uniq.gtf intersect.gtf");
system("rm noverlap.gtf overlap.gtf overlap_par.gtf overlap_par_uniq.gtf overlap_par_diff.gtf");
system("rm noverlap_count.txt overlap_count.txt chr_file1.txt chr_file2.txt snd_file1.txt");
system("rm snd_file2.txt gtf1_count.txt gtf2_count.txt chr_gtf1.txt chr_gtf2.txt chr_intern1.txt chr_intern2.txt snd_gtf1.txt snd_gtf2.txt");
system("rm overlap_trans.gtf overlap_gene.gtf overlap_exon.gtf overlap_diff1.gtf overlap_diff2.gtf");
##################################################################################################


############### User defined subroutines ######################

##############################################################################################
####This subroutine writes the hash created by 'chrcount' or 'sndcount' subroutine in
####a tab deliminated format to a text file which is then used by R-Script for plotting purposes.
##############################################################################################

sub tabout {
    my $subj = shift;
    my $prefix = shift || "";
    if ( ref($subj) eq "HASH" ){
        for ( keys %$subj ){
            tabout( $subj->{$_}, $prefix ."$_,");
        }
    }
    else{
        print  $main::f "$prefix$subj\n";
    }
    }


###################################################################################################
###This subroutine counts the number of different features in a <file> and creates a hash of it.
###################################################################################################

sub count{
	my $countfile = shift;
	open my $cf, $countfile or die "can not open file $countfile: $!";
	my $linecount =0;
	my $op_ref;
	while(my $a = <$cf>){
	chomp($a);
	next if(($a =~ m/^\#/) or (($a =~ m/^\"/))or ($a eq "")or ($a =~ /^\s*$/));
	my @row = split(/\t/, $a);
	$| = 1;
	$linecount++;
	my $chr     = $row[0];
        my $feature = $row[2];
	my $strand  = $row[6];
	$op_ref->{$feature} += 1;
	}
	return $linecount,$op_ref;
	}

###########################################################################
###This subroutine Breakdowns the counts of different features in a <file>
###according to the chromosomes in that <file> and creates the hash of it.
###########################################################################

sub chrcount{
	my $countfile = shift;
	open my $G, $countfile or die "can not open file $countfile: $!";
	my $linecount =0;
	my $op_ref;
	while(my $a = <$G>){
	chomp($a);
	next if(($a =~ m/^\#/) or (($a =~ m/^\"/))or ($a eq "")or ($a =~ /^\s*$/));
	my @row = split(/\t/, $a);
	$| = 1;
	$linecount++;
	my $chr     = $row[0];
        my $feature = $row[2];
	$op_ref->{$chr}->{$feature} += 1;
	}
	return $linecount,$op_ref;
	}

###########################################################################
###This subroutine Breakdowns the counts of different features in a <file>
###according to the strands in that <file> and creates the hash of it.
###########################################################################

sub sndcount{
	my $countfile = shift;
	open my $G, $countfile or die "can not open file $countfile: $!";
	my $linecount =0;
	my $op_ref;
	while(my $a = <$G>){
	chomp($a);
	next if(($a =~ m/^\#/) or (($a =~ m/^\"/))or ($a eq "")or ($a =~ /^\s*$/));
	my @row = split(/\t/, $a);
	$| = 1;
	$linecount++;
	my $chr     = $row[0];
        my $feature = $row[2]; # exon
	my $start   = $row[3]; #
	my $end     = $row[4]; #
	my $strand  = $row[6]; #
	$op_ref->{$strand}->{$feature} += 1;
	}
	return $linecount,$op_ref;
	}

##############################################################################################
####This subroutine writes the hash created by 'count' subroutine in a tab deliminated format
####to a text file which is then used by R-Script for plotting purposes.
##############################################################################################
sub feafile{
	my $outputfile = $_[0];
	my $filecount = $_[1];
	my $opr = $_[2];
        open my $fh, ">$outputfile" or die "can not open file $outputfile: $!";
        print $fh "total\t$filecount\n";
        while (my ($key,$value) = each %{$opr})
        {
            print $fh "$key\t$value\n";
        }
        close $fh;
	}


####################################################################################################################
#1) Check the strand for a feature in <fileA>, if it is '-' swaps the start and end position of that feature
#       if the strand for a feature in <fileB> is '-' swaps the start and end of that feature.
#2) Once all the swapping is done, find the differences between the start positions of <fileA> and <fileB>
#       and end positions of <fileA> and <fileB>.
#3) New file is created with both start and end position differences corresponsing to the feature in last 2 columns.
####################################################################################################################
sub swap_diff{
	my $ipovr = $_[0];
	my $opovr = $_[1];
        open my $ipo, $ipovr or die "can not open file $ipovr: $!";
        open my $opo, ">$opovr" or die "can not open file $opovr: $!";
        my $totline = `wc -l $ipovr | cut -d ' ' -f1`;
        my $temp=0;
        my $temp1=0;
        my $srt_diff=0;
        my $end_diff=0;
        while(my $line1 = <$ipo>){
            chomp($line1);
            next if(($line1 =~ m/^\#/) or (($line1 =~ m/^\"/))or ($line1 eq "")or ($line1 =~ /^\s*$/));
            my @rowovr = split(/\t/, $line1);
            $| = 1;
            #print "\t- Working on line $line_ovr of $totline\r" if $line_ovr % 10 == 0;
            my $A_chr = $rowovr[0];
            my $A_src = $rowovr[1];
            my $A_fea = $rowovr[2];
            my $A_srt = $rowovr[3];
            my $A_end = $rowovr[4];
	    my $A_scr = $rowovr[5];
            my $A_snd = $rowovr[6];
            my $B_chr = $rowovr[9];
            my $B_fea = $rowovr[11];
            my $B_srt = $rowovr[12];
            my $B_end = $rowovr[13];
            my $B_snd = $rowovr[15];
            my $AB_BP = $rowovr[18];

	    if($A_snd eq '-'){
                $temp=$A_srt;
                $A_srt=$A_end;
                $A_end=$temp;
                }
            if($B_snd eq '-'){
                $temp1=$B_srt;
                $B_srt=$B_end;
                $B_end=$temp1;
                }
            $srt_diff=$A_srt-$B_srt;
            $end_diff=$A_end-$B_end;
	    print $opo "\n$A_chr\t$A_src\t$A_fea\t$A_srt\t$A_end\t$A_scr\t$A_snd\t$B_chr\t$B_fea\t$B_srt\t$B_end\t$B_snd\t$AB_BP\t$srt_diff\t$end_diff";
        }
}


#############################################################################
#1) Check if first field chromosome is in proper format 'chr1...'
#       if not then make it in correct format
#2) Second field which is for Source and Ninth field which is for attributes
#       are made NULL. As it makes easier to do parsing after intersectBed
#       and these field are of no concern for comparison.
#############################################################################
sub convert{
	my $ipgtf = $_[0];
	my $opgtf = $_[1];
	open my $ip, $ipgtf or die "can not open file $ipgtf: $!";
	open my $op, ">$opgtf" or die "can not open file $opgtf: $!";#}
	while(my $line = <$ip>){
	chomp($line);
	next if(($line =~ m/^\#/) or (($line =~ m/^\"/))or ($line eq "")or ($line =~ /^\s*$/));
	my @row = split(/\t/, $line);
	$| = 1;
	my $chr;
        if ($row[0] eq "1") {$chr = "chr1";}elsif  ($row[0] eq "2") {$chr = "chr2";}elsif  ($row[0] eq "3") {$chr = "chr3";}
        elsif  ($row[0] eq "4") {$chr = "chr4";}elsif  ($row[0] eq "5") {$chr = "chr5";}elsif  ($row[0] eq "6") {$chr = "chr6";}
        elsif  ($row[0] eq "7") {$chr = "chr7";}elsif  ($row[0] eq "8") {$chr = "chr8";}elsif  ($row[0] eq "9") {$chr = "chr9";}
        elsif  ($row[0] eq "10") {$chr = "chr10";}elsif  ($row[0] eq "11") {$chr = "chr11";}elsif  ($row[0] eq "12") {$chr = "chr12";}
        elsif  ($row[0] eq "13") {$chr = "chr13";}elsif  ($row[0] eq "14") {$chr = "chr14";}elsif  ($row[0] eq "15") {$chr = "chr15";}
        elsif  ($row[0] eq "16") {$chr = "chr16";}elsif  ($row[0] eq "17") {$chr = "chr17";}elsif  ($row[0] eq "18") {$chr = "chr18";}
        elsif  ($row[0] eq "19") {$chr = "chr19";}elsif  ($row[0] eq "20") {$chr = "chr20";}elsif  ($row[0] eq "21") {$chr = "chr21";}
        elsif  ($row[0] eq "22") {$chr = "chr22";}elsif  ($row[0] eq "23") {$chr = "chr23";}elsif  ($row[0] eq "X") {$chr = "chrX";}
	elsif  ($row[0] eq "Y") {$chr = "chrY";} else {$chr = $row[0];}
        my $source      = "NULL";
	my $feature     = $row[2];
	my $start       = $row[3];
	my $end         = $row[4];
	my $score       = $row[5];
	my $strand      = $row[6];
	my $frame       = $row[7];
	my $attribute   = "NULL";
	print $op "$chr\t$source\t$feature\t$start\t$end\t$score\t$strand\t$frame\t$attribute\n";
	}
    return $opgtf;
}

##########################################################
###### Check if any file missing in input command#########
##########################################################
sub check_options {
    my $opts = shift;
    if( $opts->{'input_gtf1'} ) {
        $input_gtf1 = $opts->{'input_gtf1'};
    } else {
        die("Option input_gtf1|g1 is required.");
    }
    if( $opts->{'input_gtf2'} ) {
        $input_gtf2 = $opts->{'input_gtf2'};
    } else {
        die("Option input_gtf2|g2 is required.");
    }
}
