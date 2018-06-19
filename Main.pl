#! /usr/bin/perl
use strict;
use FindBin qw($Bin);
my $root_path = $Bin;
my $len_file = "$root_path/Database/hg19/hg19.fa.fai";

die "Usage:perl $0 <indexed sorted bam> <vcf> <window size> <max LFR> <min SupportBarcodes> <min Link> <shellDir> <outDir_tmp> <outDir_final>
Example:perl $0 input.bam input.vcf 10000000 300000 1 1 shDir outDir_tmp outDir_final

This script is writen for individual genome phasing based on hg19 reference genome.
<window size> is set for each phasing window, 10MB is suggested.
<max LFR> is decided by the maximum length of DNA fragment
<min SupportBarcodes> is minimum count of supporting barcode for each variant
<min Link> is minimum count of linked info

Author: Sun Yuhui
Email: sunyuhui\@genomics.cn\n" unless @ARGV==9;


my $bam = $ARGV[0];
my $vcf = $ARGV[1];
my $win_size = $ARGV[2];
my $max_LFR_len = $ARGV[3];
my $min_barcode_num = $ARGV[4];
my $min_link = $ARGV[5];

my $shellDir = $ARGV[6];
mkdir $shellDir unless -d $shellDir;

my $outDir = $ARGV[7];
mkdir $outDir unless -d $outDir;

my $mergeDir = $ARGV[8];
mkdir $mergeDir unless -d $mergeDir;


my $size = 10000000;
open LEN, $len_file or die $!;
while(<LEN>){
	chomp;
	my ($chr, $len) = (split /\s+/, $_)[0,1];
	next if /_/;
	next unless /chr[\d+|X]/;
	my $num = int($len/$size)+1;
	my ($beg, $end);
	open SH, ">$shellDir/run.$chr.sh" or die $!;
	for my $i (1..$num){
		if ($i eq 1){
			$beg = 1;
			$end = $size;
		}elsif ($i >1 and $i <$num){
			$beg = ($i-1)*$size+1-500000;
			$end = $i *$size;
		}elsif ($i eq $num){
			$beg = ($i-1)*$size+1-500000;
			$end = $len;
		}
		print SH "perl $root_path/Bin/LongHap.pl $bam $vcf $chr:$beg-$end $max_LFR_len $min_barcode_num $min_link $outDir/$chr\_$beg\_$end.log 2>$outDir/$chr\_$beg\_$end.error\n";
	}
	print SH "perl $root_path/Bin/merge.pl $outDir $chr $win_size $mergeDir/$chr.log\n";
	close SH;
}
close LEN;

