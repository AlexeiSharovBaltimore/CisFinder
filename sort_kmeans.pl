#!/usr/bin/perl
#Sorts the output from the kmeans_motif.c program
#sort_kmeans.pl input.txt output.txt

use strict;
my $usage = "$0 input output\n";
my $arg=0;
my $input_file = $ARGV[$arg++] or die $usage;
my $output_file = $ARGV[$arg++] or die $usage;

my %hashCl;
my @motifs;
my @clusters;
open(INFO, $input_file) or die "File $input_file not found\n";
while(my $line=<INFO>){
	$line =~ s/\s+$//;
	if(!$line){ next; }
	my @items = split(/\t/,$line);
	if($items[0] =~ /^\d+$/){
		push(@clusters,\@items);
	}elsif($items[1]){
		$hashCl{$items[1]}++;
		push(@motifs,\@items);
	}else{
		print "ERR\n";
	}
}
close INFO;

my @iClass;
my $count=1;
foreach my $key (sort {$hashCl{$b}<=>$hashCl{$a}} keys %hashCl){
	$iClass[$key] = $count++;
}
open (OUT, ">$output_file");
foreach my $ref (@motifs){
	$ref->[1] = $iClass[$ref->[1]];
}
foreach my $ref (@clusters){
	$ref->[0] = $iClass[$ref->[0]];
}
foreach my $ref (sort {$a->[1]<=>$b->[1]} @motifs){
	print OUT join("\t",@$ref)."\n";
}
foreach my $ref (sort {$a->[0]<=>$b->[0]} @clusters){
	my @data = @$ref;
	splice(@data,2,0,$data[0]);
	if(!$data[0]){ next; }
	print OUT join("\t",@data)."\n";
}
close OUT;
exit(0);

#**************************************
sub  convert {
#**************************************
my $ref = shift;
my $score = shift;

my @off;
my $n;
if($ref){ $n = @$ref; }
if(!$n || !$score){
	print "empty\n"; return();
}
my @out;
my @len;
my $global_offset=0;
my $global_length=0;
for(my $i=0; $i<$n; $i++){
	my ($score1,$off) = compare_patterns($ref->[0],$ref->[$i]);
	push(@off,$off);
	$score->[$i] = $score1;
	if($off<0 && $global_offset < -$off){ $global_offset = -$off; }
	if($global_length < $off+length($ref->[$i])){ $global_length = $off+length($ref->[$i]) }
	#print "A $ref->[$i]\t$ref->[0] $score1,$off\n";
}
for(my $i=0; $i<$n; $i++){
	my $nn = $off[$i]+$global_offset;
	my $space = "";
	for(my $j=0; $j<$nn; $j++){
		$space .= " ";
	}
	$nn = $global_length-$off[$i]-length($ref->[$i]);
	my $revPattern = reverse_pattern($ref->[$i]);
	$ref->[$i] = $space.$ref->[$i];
	$space = "";
	for(my $j=0; $j<$nn; $j++){
		$space .= " ";
	}
	$ref->[$i] .= "\t".$space.$revPattern;
}
return();
}

#**************************************
sub  reverse_pattern{
#**************************************
my $pattern = shift;
my $revPattern = reverse($pattern);
$revPattern =~ tr/ACGTRYMKacgtrymk/TGCAYRKMtgcayrkm/;
return($revPattern);
}

#**************************************
sub  compare_patterns{
#**************************************
my $p1 = shift;
my $p2 = shift;

my $n1 = length($p1);
my $n2 = length($p2);
my @pw1 = split(//,$p1);
my @pw2 = split(//,$p2);
my $offset_start = -10;
my $offset_end = $n1-$n2+10;
my $off_out=0;
if($offset_end > 10){ $offset_end=10; }
if($n2>$n1){
	$offset_start = $n1-3-$n2;
	if($offset_start < -10){ $offset_start=-10; }
	$offset_end = 8;
}
my $score=-100;
for(my $off=$offset_start; $off<$offset_end; ++$off){
	my ($start,$end,$sum)=(0,0,0);
	my ($overstart,$overend) = (0,$n1);
	if($off<0){ $start = $off; }
	else{ $overstart = $off; }
	if($off+$n2 < $n1){ $overend = $off+$n2; }
	my $overlap = $overend-$overstart;
	if($overlap < 6){ next; }
	my $n = $n1+$n2-$overlap;
	$end = $n1;
	if($off+$n2 > $n1){ $end=$off+$n2; }
	for(my $i1=$start; $i1<$end; ++$i1){
		my $i2 = $i1-$off;
		my $i3 = $n2-$i2-1;
		$n++;
		if($i1>=0 && $i1<$n1 && $i2>=0 && $i2<$n2){
			my $pw11 = $pw1[$i1];
			my $pw22 = $pw2[$i2];
			if($pw11 eq $pw22){
				if($pw11=~/[A-Z]/){ $sum += 4.0; }
				elsif($pw11=~/[a-z]/){ $sum += 3.0; }
				else{ $sum += 1; }
			}elsif(lc($pw11) eq lc($pw22)){
				if(lc($pw11)=~/[a-z]/){ $sum += 2.5; }
				else{ $sum += 1; }
			}else{
				if(lc($pw11) eq "w" && lc($pw22)=~/[at]/ ||
				   lc($pw22) eq "w" && lc($pw11)=~/[at]/ ||
				   lc($pw11) eq "s" && lc($pw22)=~/[cg]/ ||
				   lc($pw22) eq "s" && lc($pw11)=~/[cg]/){
					if($pw11=~/[A-Z]/ && $pw22=~/[A-Z]/){ $sum += 1.5; }
					else{ $sum += 0.4; }
				}
				elsif($pw11=~/[ACGT]/ && $pw22=~/[ACGT]/){ $sum -= 4; }
				elsif($pw11=~/[acgtsw]/ && $pw22=~/[acgtsw]/){ $sum -= 2; }
				elsif($pw11=~/[ACGT]/ && $pw22 ne "_" || $pw22=~/[ACGT]/ && $pw11 ne "_"){ $sum -= 2; }
				else{ $sum -= 1.5; }
			}
		}
	}
	if(!$n){ next; }
	my $score1 = $sum/$n;
	#print "C $off  $score1 $sum $n\n";  
	if($score1 > $score){
		$score = $score1;
		$off_out = $off;
	}
}
return ($score,$off_out);
}


