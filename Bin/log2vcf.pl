#! /usr/bin/perl
use strict;

my %hap1;
my %hap2;
my %block;
$/ = "break";
open LOG, $ARGV[0] or die $!;
<LOG>;
while(<LOG>){
	my ($line1, $line2) = (split /\n+/, $_)[1,2];
	my @tmp1 = split /;/, $line1;
	my @tmp2 = split /;/, $line2;
	my $num = @tmp1;
	if ($num eq 1){
		next;
	}elsif ($num >1){
		my $min;
                my $max;
		my ($chr, $pos,$alle);
                for (my $i=0; $i<@tmp1; $i++){
                        $tmp1[$i] =~ /(.+)_(\d+)_(\w+)/;
                        ($chr, $pos,$alle) = ($1, $2,$3);	
                        $min = $pos if $i==0;
                        $min = $pos if $pos <$min;
                        $max = $pos if $pos >$max;
		}	
		my $block_id = "$chr\_$min\_$max";
		for (my $i=0; $i<@tmp1; $i++){
			$tmp1[$i] =~ /(.+)_(\d+)_(\w+)/;
			($chr, $pos, $alle) = ($1, $2, $3);
			$hap1{"$chr:$pos"} = $alle;
			$block{"$chr:$pos"} = $block_id;
		}
		for (my $i=0; $i<@tmp2; $i++){
                        $tmp2[$i] =~ /(.+)_(\d+)_(\w+)/;
                        ($chr, $pos, $alle) = ($1, $2, $3);
                        $hap2{"$chr:$pos"} = $alle;
			$block{"$chr:$pos"} = $block_id;
                }	
	}
}
close LOG;
$/ = "\n";

my %het;
open HET, $ARGV[1] or die $!;
while(<HET>){
	my ($chr, $pos, $alle, $barcodes) = (split /\s+/, $_)[2,3,4,5];
	$het{"$chr:$pos"}{$alle} =$barcodes;
}
close HET;

my %hom;
open HOM, $ARGV[2] or die $!;
while(<HOM>){
	my ($chr, $pos, $alle, $barcodes) = (split /\s+/, $_)[2,3,4,5];
	$hom{"$chr:$pos"}{$alle} =$barcodes;
}
close HOM;

open VCF, $ARGV[3] or die $!;
while(<VCF>){
	chomp;
	if (/^#/){
		print $_,"\n";
		if (/^##FORMAT=<ID=PL/){ ###FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
			print "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Block ID, HETVAR stands for unphased heterozygous variants while HOMVAR stands for homozygous variants\">\n";
			print "##FORMAT=<ID=AL1,Number=1,Type=String,Description=\"Allele in Hap1\">\n";
			print "##FORMAT=<ID=AL2,Number=1,Type=String,Description=\"Allele in Hap2\">\n";
			print "##FORMAT=<ID=BA1,Number=1,Type=String,Description=\"Barcode IDs supporting AL1\">\n";
			print "##FORMAT=<ID=BA2,Number=1,Type=String,Description=\"Barcode IDs supporting AL2\">\n";
		}	
	}else{
		my @tmp = split /\s+/, $_;
		my ($chr, $pos) = @tmp[0,1];
		$pos = "$chr:$pos"; #attention
		my ($format, $detail) = @tmp[8,9]; #GT:AD:DP:GQ:PL  1/1:0,5:5:21:251,21,0
		$format .= ":BI:AL1:AL2:BA1:BA2";
		if ($het{$pos}){                 #hete
			if ($block{$pos}){ #phased
				my $alle_1 = $hap1{$pos};
				my $alle_2 = $hap2{$pos};
				my $barcodes_1 = $het{$pos}{$alle_1};
				my $barcodes_2 = $het{$pos}{$alle_2};
				$detail .= ":$block{$pos}:$alle_1:$alle_2:$barcodes_1:$barcodes_2";
			}else{                    #unphased
				$detail .= ":HETVAR::::";
			}
			for my $i (0..7){
				print $tmp[$i],"\t";
			}
			print $format, "\t", $detail, "\n";
		}elsif ($hom{$pos}){ 		#homo
			$detail .= ":HOMVAR::::";
			for my $i (0..7){
				print $tmp[$i],"\t";
            		}
            		print $format, "\t", $detail, "\n";
		}else{                      		#filtered
			next;
		}
	}
}
close VCF;
