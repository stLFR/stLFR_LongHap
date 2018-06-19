#! /usr/bin/perl
use strict;
use FindBin qw($Bin); 
my $root_path = $Bin;

die "Usage:perl $0 <indexed sorted bam> <vcf> <region> <max LFR> <min SupportBarcodes> <min Link> <output>
Example:perl $0 sorted.rmdup.bam final.vcf chr1:1-10000000 300000 3 1 chr1_1_10000000.phased

<region> is a specific region in which the variants will be phased
<max LFR> is decided by the maximum length of DNA fragment
<min SupportBarcodes> is minimum count of supporting barcode for each variant
<min Link> is minimum count of linked info

Author: Sun Yuhui
Email: sunyuhui\@genomics.cn\n" unless @ARGV==7;

my $bam = $ARGV[0];
my $var = $ARGV[1];
my $region = $ARGV[2];
my $maxLFR = $ARGV[3];
my $minSUPPORT= $ARGV[4];
my $minLink = $ARGV[5];
my $output = $ARGV[6];
my $outDir;
if ($output =~ /(.+)\/.+/){
	$outDir = $1;
}else{
	$outDir = ".";
}


$region =~ /(.+):(\d+)-(\d+)/;
my ($chr, $beg, $end) = ($1, $2, $3);

my (%marker, %marker_hete2, %marker_hete1, %marker_homo);
my %het;
my %var;
my $num_ori;
open VAR, $var or die $!;
while(<VAR>){
	chomp;
	next if /^#/;
	my ($ref, $pos, $ref_base, $alt_base) = (split /\s+/, $_)[0,1,3,4];
	
	next if $ref ne $chr;
	next unless ($pos >= $beg and $pos <= $end);
	
	next if /LowQual/; #discard low quality variants
	
#	next unless length $ref_base == 1; #part of indels, especially deletions
	my $info = (split /\s+/, $_)[-1];
	my $depth = (split /:/, $info)[2];
	if ($info=~/1\/2/){ #het-het
		$alt_base =~ /(\w+),(\w+)/;
		my ($alle_1, $alle_2) = ($1, $2);
		my $new_alle_1 = alle ($alle_1, $ref_base);
		my $new_alle_2 = alle ($alle_2, $ref_base);
		$marker{$pos} = "$new_alle_1:$new_alle_2";
		$marker_hete2{$pos} = "$new_alle_1:$new_alle_2";
		$var{$pos}{$new_alle_1}=1;
                $var{$pos}{$new_alle_2}=1;
                $het{$pos}=1;
                $num_ori++;
	}elsif ($info=~/0\/1/){ #het-ref
		my $new_alle = alle ($alt_base, $ref_base);
		$ref_base =~ /^([A|T|G|C])/;
		$ref_base = $1;
		$marker{$pos} = "$ref_base:$new_alle";
                $marker_hete1{$pos} = "$ref_base:$new_alle";
                $var{$pos}{$new_alle}=1;
                $het{$pos}=1;
                $num_ori++;
	}elsif ($info=~/1\/1/ or $info=~/1|1/){ #homo
		my $new_alle = alle ($alt_base, $ref_base);
		$ref_base =~ /^([A|T|G|C])/;
                $ref_base = $1;
		$marker{$pos} = "$ref_base:$new_alle";
                $marker_homo{$pos} = "$ref_base:$new_alle";
	}else{ #missed other type?
		die "missed other type:\n$_\n";
	}
}
close VAR;
warn "Reading VCF finished!\n";

my $min;
foreach my $pos (sort {$a<=>$b} keys %marker){
	$min= $pos;
	last;
}

my $reads_num;  my %info; my %final;
open IN, "$root_path/samtools-1.3/bin/samtools view $bam $region | " or die $!;
while (<IN>){
	chomp;
	next if /^@/;
	$reads_num++;
	if (($reads_num % 100000) == 0){
		warn "$reads_num reads finished!\n";
	}
	my ($qname,$pos,$cigar,$seq)=(split /\s+/, $_)[0,3,5,9];
	next if $cigar eq "*"; #discard unmapped read
	next if $qname =~ /#0_0_0$/; #discard reads of unknown barcode
	my $ac_len = 0;
	while ($cigar !~ /^$/){
		if ($cigar =~ /^([0-9]+[MIDSH])/){
			my $cigar_part = $1;
			if ($cigar_part =~ /(\d+)M/){ #match or mismatch
				$ac_len += $1;
			}elsif ($cigar_part =~ /(\d+)D/){
				$ac_len += $1;
			}
			$cigar =~ s/$cigar_part//;
		}else{
			die "Unexpected cigar: $qname $cigar\n";
		}
	}
	my $start = $pos;
	my $end = $pos + $ac_len;
	my $var = 0;
	if ($min < $start){
		delete $marker{$min};
		foreach my $rpos (sort {$a<=>$b} keys %marker){
			if ($rpos < $start){
				delete $marker{$rpos};
				next;
			}elsif ($rpos >= $start and $rpos <= $end){
				$min = $rpos;
				$var = 1;
				last;
			}elsif ($rpos > $end){
				$min = $rpos;
				last;
			}
		}
	}elsif ($min >= $start and $min <= $end){
		$var =1;
	}else{
		next;
	}
	
	next unless $var eq 1;
	
	$ac_len =0;
	my ($m_len, $ins_len, $del_len, $soft_len) = (0,0,0,0,0);
	$qname =~ /(.+)\#(.+)/;
	my $barcode = $2;
	$cigar = (split /\s+/, $_)[5];
	while ($cigar !~ /^$/){
		if ($cigar =~ /^([0-9]+[MIDSH])/){
			my $cigar_part = $1;
			if ($cigar_part =~ /^(\d+)M/){ #match or mismatch, and not an indel
				my $cigar_part_len = $1;
				$start = $pos +$ac_len;
				$end = $pos + $ac_len + $cigar_part_len-1;
				for my $i ($start..$end){
					if ($marker{$i}){
						next if ($i eq $end and $cigar_part !~ /^\d+M\d+[I|D]/);
						$marker{$i} =~ /(\w+):(\w+)/;
						my ($alle_1, $alle_2) = ($1, $2);
						my $qpos = $i - $pos +1 - $del_len + $ins_len + $soft_len;
						my $query_base = substr($seq, $qpos-1,1);
						if ($query_base eq $alle_1 or $query_base eq $alle_2){
							if ($info{$barcode}{$i}){
								if ($info{$barcode}{$i} eq $query_base){
									next;
								}else{
									$info{$barcode}{$i} = "conflict";
#									delete $final{$i}; #TODO
								}
							}else{
								$info{$barcode}{$i} = "$query_base";
								$final{$i} = "$alle_1:$alle_2";
							}	
						}else{  #existing a mismatch, but different from the called alter allele
							next;
						}
					}else{ #passed if there is no het marker here
						next;
					}
				}
				$m_len += $cigar_part_len;
				$ac_len += $cigar_part_len;
			}elsif ($cigar_part =~ /(\d+)I/){
				my $cigar_part_len = $1;
			    	$start = $pos + $ac_len-1;
				my $query_base = "I$cigar_part_len";
				if ($marker{$start}){
					$marker{$start} =~ /(\w+):(\w+)/;
                                	my ($alle_1, $alle_2) = ($1, $2);
					if ($query_base eq $alle_1 or $query_base eq $alle_2){
						if ($info{$barcode}{$start}){
                                        		unless ($info{$barcode}{$start} eq $query_base){
                                                        	$info{$barcode}{$start} = "conflict";
							}
						}else{
							$info{$barcode}{$start} = "$query_base";
							$final{$start} = "$alle_1:$alle_2";
						}
					}	
				}
				$ins_len += $cigar_part_len;
			}elsif ($cigar_part =~ /(\d+)D/){
				my $cigar_part_len = $1;
				$start = $pos + $ac_len-1;
				my $query_base = "D$cigar_part_len";
                                if ($marker{$start}){
                                        $marker{$start} =~ /(\w+):(\w+)/;
                                        my ($alle_1, $alle_2) = ($1, $2);
                                        if ($query_base eq $alle_1 or $query_base eq $alle_2){
                                                if ($info{$barcode}{$start}){
                                                        unless ($info{$barcode}{$start} eq $query_base){
                                                                $info{$barcode}{$start} = "conflict";
                                                        }
                                                }else{
                                                        $info{$barcode}{$start} = "$query_base";
							$final{$start} = "$alle_1:$alle_2";
                                                }
                                        }
                                }
				$del_len += $cigar_part_len;
				$ac_len += $cigar_part_len;
			}elsif ($cigar_part =~ /(\d+)S/){
				my $cigar_part_len = $1;
				$soft_len += $cigar_part_len;
			}
			$cigar =~ s/$cigar_part//;
		}else{
			die "Unexpected cigar: $qname $cigar\n";
		}
	}
	unless ($m_len + $ins_len + $soft_len == length $seq){
		die "Error in seq length in $qname!\n";
	}
}
close IN;
warn "Reading BAM finished!\n";

my %support;
my %barcodes;
my %barcodes_hom;
foreach my $barcode (sort keys %info){
	foreach my $pos (sort {$a<=>$b} keys %{$info{$barcode}}){
		my $base = $info{$barcode}{$pos};
		next unless ($base =~ /[A|T|G|C]/ or $base=~ /[I|D]\d+/);
		if ($het{$pos}){
			if ($barcodes{"$chr:$pos:$base"}){
				$barcodes{"$chr:$pos:$base"} .= ";$barcode";
			}else{
				$barcodes{"$chr:$pos:$base"} = "$barcode";
			}
			$support{$pos}++ if $var{$pos}{$base}; #add one, number of barcodes supporting this variant
		}else{
			if ($barcodes_hom{"$chr:$pos:$base"}){
				$barcodes_hom{"$chr:$pos:$base"} .= ";$barcode";
			}else{
				$barcodes_hom{"$chr:$pos:$base"} = "$barcode";
			}
		}
	}
}


open HET, ">$outDir/$chr\_$beg\_$end.hete.barcodes" or die $!;
open HOM, ">$outDir/$chr\_$beg\_$end.homo.barcodes" or die $!;
open UN, ">$outDir/$chr\_$beg\_$end.unknown.barcodes" or die $!;
foreach my $pos (keys %marker_hete2){
	$marker_hete2{$pos} =~ /(\w+):(\w+)/;
	my ($alt_1, $alt_2) = ($1, $2);
	if ($barcodes{"$chr:$pos:$alt_1"} and $barcodes{"$chr:$pos:$alt_2"}){
		print HET "NULL\tNULL\t$chr\t$pos\t$alt_1\t",$barcodes{"$chr:$pos:$alt_1"},"\n";
		print HET "NULL\tNULL\t$chr\t$pos\t$alt_2\t",$barcodes{"$chr:$pos:$alt_2"},"\n";
	}else{
                print UN "NULL\tNULL\t$chr\t$pos\t$alt_1\t$alt_2\n";
        }
}
foreach my $pos (keys %marker_hete1){
        $marker_hete1{$pos} =~ /(\w+):(\w+)/;
        my ($ref_alle, $alt_alle) = ($1, $2);
        if ($barcodes{"$chr:$pos:$alt_alle"}){
                print HET "NULL\tNULL\t$chr\t$pos\t$ref_alle\t",$barcodes{"$chr:$pos:$ref_alle"},"\n";
                print HET "NULL\tNULL\t$chr\t$pos\t$alt_alle\t",$barcodes{"$chr:$pos:$alt_alle"},"\n";
        }else{
                print UN "NULL\tNULL\t$chr\t$pos\t$ref_alle\t$alt_alle\n";
        }
}
foreach my $pos (keys %marker_homo){
        $marker_homo{$pos} =~ /(\w+):(\w+)/;
        my ($ref_alle, $alt_alle) = ($1, $2);
        if ($barcodes_hom{"$chr:$pos:$alt_alle"}){
                print HOM "NULL\tNULL\t$chr\t$pos\t$ref_alle\t",$barcodes_hom{"$chr:$pos:$ref_alle"},"\n";
                print HOM "NULL\tNULL\t$chr\t$pos\t$alt_alle\t",$barcodes_hom{"$chr:$pos:$alt_alle"},"\n";
        }else{
                print UN "NULL\tNULL\t$chr\t$pos\t$ref_alle\t$alt_alle\n";
        }
}
close HET;
close HOM;
close UN;

my %link;
foreach my $barcode (sort keys %info){
	my %tmp;
	foreach my $pos (sort {$a<=>$b} keys %{$info{$barcode}}){
		my $base = $info{$barcode}{$pos};
		next unless ($base =~ /[A|T|G|C]/ or $base=~ /[I|D]\d+/);
		next unless $support{$pos} >= $minSUPPORT;
		foreach my $fpos (keys %tmp){
			my $fbase = $tmp{$fpos};
			if ((abs ($pos - $fpos)) < $maxLFR){
				$link{$pos}{$fpos}{"$base-$fbase"} .= "$barcode;";
				$link{$fpos}{$pos}{"$fbase-$base"} .= "$barcode;";
			}
		}
		$tmp{$pos}=$base;
	}
}

open LINK, ">$outDir/$chr\_$beg\_$end.link.list" or die $!;
foreach my $pos (sort {$a<=>$b} keys %link){
	foreach my $rpos (sort {$a<=>$b} keys %{$link{$pos}}){
		foreach my $trans (keys %{$link{$pos}{$rpos}}){
			print LINK "$pos\t$rpos\t$trans\t$link{$pos}{$rpos}{$trans}\n";
		}
	}	
}
close LINK;

my ($seed_1, $seed_2, %add, %danger);
foreach my $pos (sort {$a<=>$b} keys %final){
        next unless $support{$pos} >= $minSUPPORT;
        $final{$pos} =~ /(\w+):(\w+)/;
        my ($geno_1, $geno_2) = ($1, $2);
        $seed_1 = "$chr\_$pos\_$geno_1";
        $seed_2 = "$chr\_$pos\_$geno_2";
        $add{$pos} =1;
        last;
}

open OUT, ">$output" or die $!;
open LEN, ">$outDir/$chr\_$beg\_$end.len.list" or die $!;
open PHASED, ">$outDir/$chr\_$beg\_$end.phased.var.barcode" or die $!;
open UNPHASED, ">$outDir/$chr\_$beg\_$end.unphased.var.barcode" or die $!;
open CONFLICT, ">$outDir/$chr\_$beg\_$end.conflict.list" or die $!;
my $num_block;
my $num_phased_marker;
my (%merge1, %merge2);
my $max;
while($seed_1){
        my $extend =0;
        foreach my $pos (sort {$a<=>$b} keys %link){
                last if ($pos-$max>500000);
                next if $add{$pos};
                $final{$pos} =~ /(\w+):(\w+)/;
                my ($geno_a, $geno_b) = ($1, $2);
                my ($yes_p1, $yes_p2, $no) = (0,0,0);
                foreach my $locus_pos (keys %{$link{$pos}}){
                        next unless $merge1{$locus_pos};
                        my ($geno_1, $geno_2) = ($merge1{$locus_pos}, $merge2{$locus_pos});
#                	next unless ((length $geno_1) ==1 and (length $geno_2) ==1); #ignore indel within the seeds when phasing
		        my ($a_1, $b_2, $a_2, $b_1) = (0,0,0,0);
                        foreach my $trans (keys %{$link{$locus_pos}{$pos}}){
                                my $num = $link{$locus_pos}{$pos}{$trans};
                                $trans =~ /(\w+)-(\w+)/;
				if ($geno_1 eq $1){
                                	if ($geno_a eq $2){
                                               	$a_1 += $num;
	                                }elsif ($geno_b eq $2){
        	                                $b_1 += $num;
					}
                        	}elsif ($geno_2 eq $1){
                                	if ($geno_b eq $2){
                                  		$b_2 += $num;
	                                }elsif ($geno_a eq $2){
        	                        	$a_2 += $num;
                	        	}
                        	}
                        }
                        my $p1 = $a_1+ $b_2;
                        my $p2 = $a_2+ $b_1;
                        if ( $p1 >= $minLink and $p1 > 5*$p2){ #TODO
                                $yes_p1++;
                        }elsif ($p2 >= $minLink and $p2 > 5*$p1){ #TODO
                                $yes_p2++;
                        }elsif($p1 >= $minLink and $p2 >= $minLink and ($p1/($p1+$p2)) >0.3 and ($p1/($p1+$p2)) <0.7){ #TODO
                                $no++;
				print CONFLICT $chr, "\t", $pos, "\n";
                        }
                }
                my ($seed_1_geno, $seed_2_geno);
                if ($yes_p1 eq $yes_p2){
                        print CONFLICT $chr, "\t", $pos, "\n";
			next;
                }else{
                        if ($yes_p1>$yes_p2){
                                $seed_1_geno = $geno_a;
                                $seed_2_geno = $geno_b;
                        }else{
                                $seed_1_geno = $geno_b;
                                $seed_2_geno = $geno_a;
                        }
                        $seed_1 .= ";$chr\_$pos\_$seed_1_geno";
                        $seed_2 .= ";$chr\_$pos\_$seed_2_geno";
                        $add{$pos}=1;
                        $max = $pos if $pos >$max;
                        $num_phased_marker++;
                        $extend = 1;
                        $merge1{$pos}=$seed_1_geno;
                        $merge2{$pos}=$seed_2_geno;
                        last;
                }
        }
        if ($extend ==0){
        	print OUT "break\n",$seed_1, "\n", $seed_2, "\n";
                my @loci_1 = split /;/, $seed_1;
                my @loci_2 = split /;/, $seed_2;
                my ($seed_1_reads, $seed_2_reads);
                my $min_locus_pos=0;
                my $max_locus_pos=0;
                my $var_num;
                for (my $i=0; $i<@loci_1; $i++){
                        $var_num++;
                        my $locus_1 = $loci_1[$i];
                        $locus_1 =~ /\w+_(\w+)_(\w+)/;
                        my $locus_pos = $1;
                        $max_locus_pos = $locus_pos if $locus_pos >$max_locus_pos;
                        $min_locus_pos = $locus_pos if ($locus_pos <$min_locus_pos or $min_locus_pos==0);
                }
                my $len=$max_locus_pos-$min_locus_pos;
                print LEN $chr, "\t", $min_locus_pos, "\t", $max_locus_pos, "\t",$len,"\t", $var_num,"\n";
                if ($len > 0){
                        $num_phased_marker++; #add the start of this seed
                        $num_block++;
                        for (my $i=0; $i<@loci_1; $i++){
                                my $locus_1 = $loci_1[$i];
                                my $locus_2 = $loci_2[$i];
                                $locus_1 =~ /\w+_(\w+)_(\w+)/;
                                my ($pos, $base_1) = ($1,$2);
                                $locus_2 =~ /\w+_(\w+)_(\w+)/;
                                my $base_2 = $2;
                                print PHASED "block$num_block\thap1\t$chr\t$pos\t$base_1\t",$barcodes{"$chr:$pos:$base_1"},"\n";
                                print PHASED "block$num_block\thap2\t$chr\t$pos\t$base_2\t",$barcodes{"$chr:$pos:$base_2"},"\n";
                        }
                }else{
                        $loci_1[0] =~ /\w+_(\w+)_(\w+)/;
                        my ($pos, $base_1) = ($1,$2);
                        $loci_2[0] =~ /\w+_(\w+)_(\w+)/;
                        my $base_2 = $2;
                        print UNPHASED "NULL\tNULL\t$chr\t$pos\t$base_1\t",$barcodes{"$chr:$pos:$base_1"},"\n";
                        print UNPHASED "NULL\tNULL\t$chr\t$pos\t$base_2\t",$barcodes{"$chr:$pos:$base_2"},"\n";
                }

                undef $seed_1; #TODO
                undef $seed_2; #TODO
                undef %merge1;
                undef %merge2;
                foreach my $pos (sort {$a<=>$b} keys %final){
                        next if $add{$pos};
                        next if $danger{$pos};
                        next unless $support{$pos}>=$minSUPPORT;
                        $final{$pos}=~ /(\w+):(\w+)/;
                        my ($geno_1,$geno_2) = ($1, $2);
			$seed_1 = "$chr\_$pos\_$geno_1";
                        $seed_2 = "$chr\_$pos\_$geno_2";
                        $add{$pos}=1;
                        $max = $pos;
                        $merge1{$pos} = $geno_1;
                        $merge2{$pos} = $geno_2;
                        last;
		}
	}
}
close LEN;
close PHASED;
close UNPHASED;
close CONFLICT;
close OUT;

sub alle {
	my ($alt, $ref) = @_;
	my $alt_len = length $alt;
	my $ref_len = length $ref;
	if ($alt_len < $ref_len){
		my $indel_len = $ref_len - $alt_len;
                $alt = "D$indel_len";
	}elsif ($alt_len > $ref_len){
		my $indel_len = $alt_len - $ref_len;
		$alt = "I$indel_len";
	}else{
        	$alt =~ /^([A|T|G|C])/; #GA G,AA
		$alt = $1;
	}
	return $alt;
}

