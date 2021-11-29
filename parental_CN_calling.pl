#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';
use List::Util qw( min max );

#INPUT

my @line = ();
my @chromnum = ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36');
my @chroms = ('LmjF.01','LmjF.02','LmjF.03','LmjF.04','LmjF.05','LmjF.06','LmjF.07','LmjF.08','LmjF.09','LmjF.10','LmjF.11','LmjF.12','LmjF.13','LmjF.14','LmjF.15','LmjF.16','LmjF.17','LmjF.18','LmjF.19','LmjF.20','LmjF.21','LmjF.22','LmjF.23','LmjF.24','LmjF.25','LmjF.26','LmjF.27','LmjF.28','LmjF.29','LmjF.30','LmjF.31','LmjF.32','LmjF.33','LmjF.34','LmjF.35','LmjF.36');
my @samples = ('S1','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19','S2','S20','S21','S22','S23','S24','S25','S26','S27','S28','S29','S3','S30','S31','S32','S33','S34','S35','S36','S37','S38','S39','S4','S40','S41','S42','S5','S6','S7','S8','S9');
my $parenta = 'S1';
my $parentb = 'S3';
my $parentacol = '9';
my $parentbcol = '31';
#S1xS2 F1 samples
#my @fones = ('9','20','15','16','17','26','27','28','29','42','44','45','46','47','48','49');
#S1xS3 F1 samples
my @fones = ('9','31','10','11','12','13','14','19','21','22','35','36','37','38','39','40','41','43','50');
#S2xS3 F1 samples
#my @fones = ('20','31','18','23','24','25','30','32','33','34');

my $file='';

my $a=0;

#$file = 'parental_freqs/chr'.($a+1).'/'.$samples[$num];
#open J, ">$file";

#for ($a=0; $a<@chroms; $a++) {
my $vcf='combined.relaxedFilter.biallelic.noMissing.vcf';
my $ploidyfile='concatPloidyMatrix.txt';
	
my $c=0;
for ($c=0; $c<@fones; $c++) {
	my $num = '';
	my $gta = '';
	my $gtb = '';
	my $chr = '';
	my $pos = '';
	my $col = '';
	my $dpa = '';
	my $dpb = '';
	
	$num = ($fones[$c] - 9);
	$col = $fones[$c];

	$file = $samples[$num] . '_' . $parenta . 'x' . $parentb . '_parental_CN.txt';
	open J, ">$file";
	print J "Chr\tWinStart\tWinEnd\t$parenta\t$parentb\n";

	open (H, $vcf);
	while (<H>) {
		chomp;
  		last if m/^$/;
		@line = split;
		next if ($line[0] =~ m'#');
		$gta = ((split ':', $line[$parentacol])[0]);
		$gtb = ((split ':', $line[$parentbcol])[0]);
		$dpa = ((split ':', $line[$parentacol])[2]);
		$dpb = ((split ':', $line[$parentbcol])[2]);
		$chr = $line[0];
		$pos = $line[1];
		
		my $ad = '';
		my $ada = '';
		my $adb = '';
		my $afa = '';
		my $afb = '';
		
		if (($dpa >= 20) && ($dpb >=20)) {
			if ((($gta eq '0/0') && ($gtb eq '1/1')) || (($gta eq '1/1') && ($gtb eq '0/0'))) {
				if (($gta eq '0/0') && ($gtb eq '1/1')) {
					$ad = ((split ':', $line[$col])[1]);
					$ada = ((split ',', $ad)[0]);
					$adb = ((split ',', $ad)[1]);
				}
				elsif (($gta eq '1/1') && ($gtb eq '0/0')) {
					$ad = ((split ':', $line[$col])[1]);
					$ada = ((split ',', $ad)[1]);
					$adb = ((split ',', $ad)[0]);
				}
	
				$afa = ($ada/($ada+$adb));
				$afb = ($adb/($ada+$adb));
			
				my $ploidycol = '';
				my $ploidy = '';
			
				open (I, $ploidyfile);
				my @lineb = ();
			
				while (<I>) {
					chomp;
	  				last if m/^$/;
					@lineb = split;
					if ($lineb[0] eq 'SampleID') {
						my $d=0;
						for ($d=0; $d<@lineb; $d++) {
							if ($chr eq $lineb[$d]) {
								$ploidycol = $d;
							}
						}
					}
					else {
						if ($lineb[0] eq $samples[$num]) {
							$ploidy = $lineb[$ploidycol];
						}
					}
				}

				close I;
				my $cna = ($afa * $ploidy);
				my $cnb = (-1*($afb * $ploidy));
				my $winstop=$pos+1;
				my $id = $chr . '_' . $pos;
				print J "$id\t$chr\t$pos\t$winstop\t$cna\t$cnb\n";
			}
		}
	}
	close H;
	close J;
}
