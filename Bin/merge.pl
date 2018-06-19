#! /usr/bin/perl
use strict;

die "Usage:perl $0 <phasedDir> <target chromosome> <bin size> <final phased output>
Example:perl $0 split_out chr1 10000000 chr1.phased
" unless @ARGV==4;

my $chr = $ARGV[1];
my $size = $ARGV[2];
my $output = $ARGV[3];

open OUT, ">$output" or die $!;

my @files = `find $ARGV[0] -name '*.log'`;
chomp @files;


my %hash;
my $total;
foreach my $file (@files){
	my $basename = (split /\//, $file)[-1];
	$basename =~ /(.+)_(\d+)_(\d+)/;
	next unless $chr eq $1;
	$total++;	
	my $num = int(($2+500000)/$size)+1;
	open IN, $file or die $!;
	$/ = "break";
	<IN>;
	while(<IN>){
		chomp;
	        my ($line1, $line2) = (split /\n/,$_)[1,2];
        #        if ((length $line1) ne (length $line2)){
	#		die "not same lenght at\n$file\n$line1\n$line2\n";
	#	}else{
	#		unless ($line1 =~ /[A|T|G|C]$/ and $line2 =~ /[A|T|G|C]$/){
	#			die "wrong ending at $file\n$line1\n$line2\n";
	#		}
	#	}
		$hash{$num} .= "break\n$line1\n$line2\n";
	}
	$/ = "\n";
}

my $passedBlock;
my $toMergeBlock; 
$toMergeBlock = "$hash{1}" if exists $hash{1};
foreach my $i (2..$total){
	unless ($hash{$i}){ #this file is empty
		$passedBlock .= "$toMergeBlock";
		$toMergeBlock ="";
		next;
	}
	my $to_add = $hash{$i};
	unless ($toMergeBlock =~ /break/){ #last file is empty
		$toMergeBlock = "$to_add";
		next;
	}
	my $beg = ($i-1)*$size+1-500000;
	my $count;
	my (%Blocks1, %Blocks2);
        my @blocks = split /break/, $toMergeBlock;
	for (my $ii=1; $ii< @blocks; $ii++){
		my $block = $blocks[$ii];
		my ($line1, $line2) = (split /\n/, $block)[1,2];
		my @sites1 = split /;/, $line1;
		my $max;
		for (my $j=0; $j <@sites1; $j++){
			$sites1[$j] =~ /(.+)_(\w+)_(\w+)/;
			my $pos = $2;
			$max = $pos if $pos >$max;
		}	
		if ($max <$beg){
			$passedBlock .= "break$block";
		}else{
			$count++;
			$Blocks1{$count} = $line1;
			$Blocks2{$count} = $line2;
		}
	}
	unless ($count>0){ #neither of this or last file is empty, but this file is NOT overlapped with next one 
		$toMergeBlock = "$to_add";
		next;
	}
	$toMergeBlock = "";
	my $last_end = ($i-1)*$size;
	my @to_add_blocks = split /break/, $to_add;
	for (my $iii=1; $iii<@to_add_blocks; $iii++){
		my $block =$to_add_blocks[$iii];
		my ($line1, $line2) = (split /\n/, $block)[1,2];
		my @sites1 = split /;/, $line1;
		my $min;
		for (my $j=0; $j <@sites1; $j++){
			$min = $sites1[$j] if $j eq 0;
			$sites1[$j] =~ /(.+)_(\w+)_(\w+)/;
                        my $pos =  $2;
                        $min = $pos if $pos <$min;
		}
		if ($min > $last_end){
			$toMergeBlock .= "break$block";
		}else{	
			$count++;
			$Blocks1{$count} = $line1;
			$Blocks2{$count} = $line2;
		}
	}
	my $seed1 = $Blocks1{1};
	my $seed2 = $Blocks2{1};
	my %inclu;
	$inclu{1} =1;
	while ($seed1){
		my %posb; my %new;
		foreach my $key (keys %Blocks1){
			next if $inclu{$key};
			my $add1 = $Blocks1{$key};
			my $add2 = $Blocks2{$key};
			my ($p1, $p2, $add1_new, $add2_new) = cal_posb($seed1, $seed2, $add1, $add2);
			$posb{$key}{1} = $p1;
			$posb{$key}{2} = $p2;
			$new{$key}{1} = $add1_new;
			$new{$key}{2} = $add2_new;
		}
		my $max=0; my $max_block;
		foreach my $key (keys %posb){
			my $overlap_num = $posb{$key}{1} + $posb{$key}{2};
			if ($overlap_num > $max){
				$max = $overlap_num;
				$max_block = $key;
			}
		}
		if ($max >0){
			my ($seed1_add, $seed2_add);
			if ($posb{$max_block}{1} >= $posb{$max_block}{2}){
				$seed1_add = $new{$max_block}{1}; 
				$seed2_add = $new{$max_block}{2};
			}else{
				$seed1_add = $new{$max_block}{2};
				$seed2_add = $new{$max_block}{1};
			}
			$seed1 .= "$seed1_add";
			$seed2 .= "$seed2_add";
			$inclu{$max_block} =1;
		}else{
			$toMergeBlock .= "break\n$seed1\n$seed2\n";
			undef $seed1;
                	undef $seed2;
			foreach my $key (keys %Blocks1){
				next if $inclu{$key};
				$seed1 = $Blocks1{$key};
				$seed2 = $Blocks2{$key};
				$inclu{$key} =1; 
				last;
			}
		}		
	}
}
print OUT $passedBlock, $toMergeBlock;
close OUT;

sub cal_posb {
	my ($seed1, $seed2, $add1, $add2) = @_;
	my ($p1, $p2) = (0,0);
	my ($add1_new, $add2_new);
	my @line1 = split /;/, $seed1;
	my @line2 = split /;/, $seed2;
	my @add_line1 = split /;/, $add1;
	my @add_line2 = split /;/, $add2;
	my %info1; my %info2;
	for (my $i=0; $i<@line1; $i++){
		$line1[$i] =~ /(.+)_(\d+)_(\w+)/;
		my ($pos, $alle_1) = ($2, $3);
		$line2[$i] =~ /(.+)_(\d+)_(\w+)/;
		my $alle_2 = $3;
		$info1{$pos} = $alle_1;
		$info2{$pos} = $alle_2;
	}
	for (my $i=0; $i<@add_line1; $i++){
		$add_line1[$i] =~ /(.+)_(\d+)_(\w+)/;
		my ($pos, $alle_1) = ($2, $3);
		$add_line2[$i]  =~ /(.+)_(\d+)_(\w+)/;
                my $alle_2 = $3;
#		if ((length $alle_1)==1 and (length $alle_2)==1){
			if ($info1{$pos} eq $alle_1 and $info2{$pos} eq $alle_2){
				$p1++;
			}elsif ($info1{$pos} eq $alle_2 and $info2{$pos} eq $alle_1){
				$p2++;
			}else{
				$add1_new .= ";$add_line1[$i]";
                        	$add2_new .= ";$add_line2[$i]";
                        	next;
			}
#		}else{
#			$add1_new .= ";$add_line1[$i]";
#			$add2_new .= ";$add_line2[$i]";
#			next;
#		}
	}
	return ($p1, $p2, $add1_new, $add2_new);
}
