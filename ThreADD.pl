#!/usr/bin/perl -w
use strict;
use FileHandle;
use List::Util qw(first);
use List::Util qw(min max);
use List::Util qw(shuffle);
use Chart::Gnuplot;


# Morag Lewis
# morag.lewis@kcl.ac.uk

#Threshold Average Difference Detection
#Perl script which compares audiometric threshold averages between samples grouped by genotype from a variant input file


##############################################################################
## Copyright 2020 Morag Lewis                                               ##
##                                                                          ##
## Licensed under the Apache License, Version 2.0 (the "License");          ##
## you may not use this file except in compliance with the License.         ##
## You may obtain a copy of the License at                                  ##
##                                                                          ##
##   http://www.apache.org/licenses/LICENSE-2.0                             ##
##                                                                          ##
## Unless required by applicable law or agreed to in writing, software      ##
## distributed under the License is distributed on an "AS IS" BASIS,        ##
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. ##
## See the License for the specific language governing permissions and      ##
## limitations under the License.                                           ##
##                                                                          ##
##############################################################################

#This script takes as input:
#1. A tab-separated file containing thresholds in the following format:
#   A header line starting "SEX" and containing sex, person ID, secondary ID, collection date, age at collection, noise history, right ear thresholds (for 8 pure tones in Hz), space, left ear thresholds (for the same 8 pure tones)
#   Values for the SEX column can be anything, but only "M" and "F" values will be used when assigning individuals to sex-specific groups. All IDs which are neither "M" nor "F" will be assigned to the "N" group, which is not used for sex-specific comparison.
#   Any unavailable threshold values should be "NA". Any other unavailable values should be ".". The script only uses sex, age and ID in addition to the tone thresholds.
#   The default 8 tones are 0.25, 0.5, 1, 2, 3, 4, 6 and 8kHz. If the first tone is 125Hz, then the script switches to use 0.125, 0.25, 0.5, 1, 2, 4, 6 and 8kHz. No other options are provided.
#2. A vcf file with variants. The IDs in the vcf header should correspond to the IDs in the threshold file, but it is not a problem if there are extra IDs in either file. Note that this script does not handle multiallelic sites, but vcf files can be "flattened" like so:
#   G 	A,C 	0/1 	0/2 	1/2
#   becomes:
#   G 	A 	0/1 	0/v 	1/v
#   G 	C 	0/v 	0/1 	v/1
#   "v" here means any allele which is neither the reference nor the alternate
#   "Minority calls" (referred to in comments) are calls which are not simple ones, eg (0/1, 1/1), and indicate a failed call
#   The script will handle mitochondrial variants, but expects them to be marked "HET1", "HOM1", "WT" or "FAIL" (for heteroplasmic, homoplasmic, wildtype or fail)
#   If variants are plotted, a gene name from the VCF file will be used as a title for the plot. A VEP-annotated vcf file will work, or you can provide your own names in the INFO field. The script expects values to be separated by "|", and the gene name and Ensembl ID should be in the fourth and fifth entry.
#   For example:
#   1       1485777 rs1622213       G       A       100.5 PASS    CSQ=A|splice_region_variant|LOW|ATAD3B|ENSG00000160072
#3. A desired minimum threshold difference, in dB
#4. A maximum permitted standard deviation. This is adjusted by the script to allow for more variability at high frequencies (+10) and less at low frequencies (-5).
#5. (Optional) A gene name (eg ATAD3B) or "gengr". In the case of a gene name being input, the script will generate graphs for every variant in that gene in the file. 
#   In the case of "gengr", the script will generate graphs for all variants in the file. No comparisons will be carried out. 
#   In both cases, the set threshold and standard deviation still need to be entered, but will not be used.
#
#The output is a list of variants which pass the required settings, and can be used as input for ThreADD_permute.pl. If a gene name or "gengr" is supplied as the fifth argument, image files with the graphs in will also be output.
#Please note that ImageMagick is required for the image generation.
#If graphs are generated without testing thresholds (eg by inputting a gene name or "gengr" as the fifth argument), the output text file instead contains details of each variant for which a graph was generated, its filename and the participant IDs in each grouping ("0_1_F" means female participants with a 0/1 genotype, ie heterozygotes)


($#ARGV > -1) or die("Usage: ThreADD.pl threshold-file call-file difference standard-deviation (gene/gengr)\n");
my ( $thresholds, $calls, $set, $seedstdv, $gene ) = @ARGV;

#x-axis settings, assuming the standard tones: 0.25, 0.5, 1, 2, 3, 4, 6, 8:
my @x = (1,2,3,4,4.56,5,5.56,6); #array for plotting
my $maxx = "6.5,0"; #coordinates for the line at y=0
my $xmx = 6.5; #maximum x axis value
my $start125 = 0; #binary to catch if starting at 0.125

my $generategraphs=0; #default=0 - do not generate graphs, just output list of passing variants

my $lookingforgene = 0; #binary to flag if variants from a specific gene are required
if ( $#ARGV > 3 ) { #if there is anything after the "set threshold", it's a gene name or a command to generate graphs
	$lookingforgene = 1;
	$generategraphs=1; #by definition, if looking for a gene, graphs are wanted
	chomp($gene);
	if ( $gene =~ /gengr/ ) { #if this is just the command to generate graphs, then all the variants are wanted
		$lookingforgene = 0;
	}
}

#binary - assume these are not mitochondrial variants. This is just for getting the legend right (ie using homoplasmic and heteroplasmic instead of homozygous and heterozygous)
my $ismit=0;

chomp($set);
chomp($seedstdv);

#populating the setstdv array
my @setstdv;
$setstdv[ 0 ] = $seedstdv - 5;
$setstdv[ 1 ] = $seedstdv - 5;
$setstdv[ 2 ] = $seedstdv;
$setstdv[ 3 ] = $seedstdv;
$setstdv[ 4 ] = $seedstdv + 5;
$setstdv[ 5 ] = $seedstdv + 5;
$setstdv[ 6 ] = $seedstdv + 10;
$setstdv[ 7 ] = $seedstdv + 10;
$setstdv[ 8 ] = $seedstdv - 5;
$setstdv[ 9 ] = $seedstdv - 5;
$setstdv[ 10 ] = $seedstdv;
$setstdv[ 11 ] = $seedstdv;
$setstdv[ 12 ] = $seedstdv + 5;
$setstdv[ 13 ] = $seedstdv + 5;
$setstdv[ 14 ] = $seedstdv + 10;
$setstdv[ 15 ] = $seedstdv + 10;

print STDERR "Looking for differences in thresholds greater than $set dB for variants in $calls\n";


#hash for holding sex, keyed on sample ID
my %sex;
#hash for holding age, keyed on sample ID
my %age;

#hash for holding thresholds, keyed on sample ID
my %thresholds;

#holds original thresholds without multiplication, used for plotting individual audiograms at the end
my %originalthrsh; 

#array for holding sample IDs
my @samples;

#Read in thresholds
my $fh = new FileHandle('<' . $thresholds );
while(<$fh>) {
	if ($_ =~ /^SEX/) { #header line
		my @fields = split(/\t/);
		if ( $fields[ 6 ] == 125 ) { #catch the case where the tones start at 0.125 - assume they are then 0.125, 0.25, 0.5, 1, 2, 4, 6, 8:
			$setstdv[ 2 ] = $seedstdv - 5;
			$setstdv[ 4 ] = $seedstdv;
			$maxx = "7.5, 0";
			$xmx = 7.5;
			@x = (1,2,3,4,5,6,6.56,7);
			$start125 = 1;
		}
        } elsif ($_ !~ /^\s/) { #ignore empty lines
	        my (@fields) = split(/\s+/);
		my $id = $fields[1];
		my @data = splice(@fields, 6,16);
		$originalthrsh{ $id } = [@data]; #makes new array reference
		my $boo;
		for $boo ( 0 .. $#data ) {
			if ( $data[ $boo ] ne "NA" ) {
				$data[ $boo ] = $data[ $boo ] * 10; #avoid floating point nonsense by multiplying by 10
			}
		}
		$thresholds{ $id } = [@data]; #makes new array reference
		my $years = $fields[4] * 100; #see above re floating point nonsense
		$age{ $id } = $years; #store age
		$sex{ $id } = $fields[0]; #store sex
		if ( $sex{ $id } !~ /M|F/ ) {
			$sex{ $id } = "N"; #tidy up any non-M/F sex designation
		}
	}
}
$fh->close();

#tracker
my $chr="0";

#hash for holding colours for graphs
my %colours = (
        "0_0" => "#000000", #black
        "0_0_M" => "#808080", #grey
        "0_0_F" => "#C0C0C0", #pale grey
        "0_0_N" => "#555555", #dark grey
        "0_1" => "#117733", #green
        "0_1_M" => "#44aa99", #teal
        "0_1_F" => "#999933", #olive
        "0_1_N" => "#88ccee", #cyan
        "v_1" => "#666633", #darker yellow
        "v_1_M" => "#eecc66", #light yellow
        "v_1_F" => "#ddcc77", #sand
        "v_1_N" => "#997700", #dark yellow
        "1_1" => "#332288", #indigo
        "1_1_M" => "#882255", #wine
        "1_1_F" => "#cc6677", #rose
        "1_1_N" => "#aa4499" #purple
);
my %gtypes = (
	"0_0" => "ref (all)",
	"0_0_M" => "ref (male)",
	"0_0_F" => "ref (female)",
	"MT0_0" => "reference (all)",
	"MT0_0_M" => "reference (male)",
	"MT0_0_F" => "reference (female)",
	"0_v" => "ref/oth (all)",
	"0_v_M" => "ref/oth (male)",
	"0_v_F" => "ref/oth (female)",
	"0_1" => "ref/alt (all)",
	"0_1_M" => "ref/alt (male)",
	"0_1_F" => "ref/alt (female)",
	"MT0_1" => "heteroplasmic (all)",
	"MT0_1_M" => "heteroplasmic (male)",
	"MT0_1_F" => "heteroplasmic (female)",
	"v_1" => "alt/oth (all)",
	"v_1_M" => "alt/oth (male)",
	"v_1_F" => "alt/oth (female)",
	"1_1" => "alt (all)",
	"1_1_M" => "alt (male)",
	"1_1_F" => "alt (female)",
	"MT1_1" => "homoplasmic (all)",
	"MT1_1_M" => "homoplasmic (male)",
	"MT1_1_F" => "homoplasmic (female)"
);
my %points = (
	"0_0" => "fill-circle",
	"0_0_M" => "fill-diamond",
	"0_0_F" => "fill-triangle",
	"0_0_N" => "fill-square",
	"0_v" => "fill-circle",
	"0_v_M" => "fill-diamond",
	"0_v_F" => "fill-triangle",
	"0_v_N" => "fill-square",
	"0_1" => "fill-circle",
	"0_1_M" => "fill-diamond",
	"0_1_F" => "fill-triangle",
	"0_1_N" => "fill-square",
	"v_1" => "fill-circle",
	"v_1_M" => "fill-diamond",
	"v_1_F" => "fill-triangle",
	"v_1_N" => "fill-square",
	"1_1" => "fill-circle",
	"1_1_M" => "fill-diamond",
	"1_1_F" => "fill-triangle",
	"1_1_N" => "fill-square"
);

#Order of plotting:
my @ordering = ("0_v_M", "0_v_F", "v_1_M", "v_1_F", "0_0_M", "0_0_F", "0_1_M", "0_1_F", "1_1_M", "1_1_F", "0_v", "v_1", "0_0", "0_1", "1_1");

#Read in data and add totals
my $fhc = new FileHandle('<' . $calls );
while(<$fhc>) {
        if (/^##/) { next; } #ignore header lines
	if (/^#/) { #get sample IDs (and print header)
		if ( $generategraphs == 1 ) {
			#header for text file
			print "Coordinate\tVariant ID\tReference\tAlternate\tFilename\n";
		} else { #if not generating graphs, just print the vcf header
			print;
		}
		my (@data) = split(/\s+/);
		my $bobby;
		for $bobby ( 9 .. $#data ) {
			my $name = $data[ $bobby ];
			chomp( $name );
			push( @samples, $name); #add the ID to the samples array
		}
	} else {
		my (@data) = split(/\s+/);
		my $outline = $_; #for outputting

		#chromosome tracking
		if ( $chr ne $data[0] ) { 
			$chr = $data[0];
			print STDERR "Chromosome $chr\n";
		}
		if ( $chr =~ /M/ ) { #detect mitochondrial data
			$ismit=1;
		} else {
			$ismit=0;
		}

		my $varinfo = $data[0].":".$data[1]."\t".$data[2]."\t".$data[3]."\t".$data[4];

		#get gene name for graphing (use Ensembl ID if no gene name given)
		my $tosplit = $data[7];
                chomp( $tosplit );
                $tosplit = $tosplit."|end"; #force split to work properly; 
		my @annots = split(/\|/,$tosplit);
		my $genename;
		if ( $annots[3] !~ /^$/ ) { 
			$genename = $annots[3];
		} else {
			$genename = $annots[4];
		}

		if ( ( $lookingforgene == 1) && ( $gene ne $genename ) ) { next; } #skip this variant if it is not in the gene wanted (and lookingforgene is set)

		#hashes for holding totals
		my %totals; #holds total counts, keyed on genotype and sex
		my %indivthrs; #holds thresholds for averaging at end, keyed on genotype and sex
		my %averages; #holds average values, keyed on genotype and sex
		my %avgage; #holds average age, keyed on genotype and sex
		my %stdevs; #holds standard deviation values, keyed on genotype and sex
		my %notables; #hash to hold details of groups which pass the criteria (if not looking at an individual gene)
		my %grpnms; #hash to hold names keyed by their genotype and sex

		my $bill;
		for $bill ( 9 .. $#data) {
			my $name = $samples[ ($bill-9) ]; #get sample ID
			if (!( exists($thresholds{ $name }) )) { next; } #if there are no thresholds for this ID, go to the next individual

			my @thrsh = @{$thresholds{ $name }}; #get the thresholds

			my $genotype = $data[ $bill ]; #get genotype
			$genotype =~ s/\//_/; #avoid "/" in entries
			chomp($genotype);
			if ( $genotype =~ /\./ ) { next; } #skip missed calls
			if ( $genotype =~ /\,/ ) { next; } #skip minority calls
			if ( $genotype eq "1_v" ) { $genotype = "v_1"; } #since v/1 and 1/v are the same

			#MT-specific settings - force homoplasmic to hom, heteroplasmic to het
			if ( $genotype eq "WT" ) { $genotype = "0_0"; }
			if ( $genotype eq "FAIL" ) { next; } #skip failures
			if ( $genotype eq "HET1" ) { $genotype = "0_1"; }
			if ( $genotype eq "HOM1" ) { $genotype = "1_1"; }

			#add sex to genotype if it contains "M" or "F"
			my $genotypex = $genotype."_N";
			if ( $sex{ $name } =~ /M|F/ ) {
				$genotypex = $genotype."_".$sex{ $name };
			}

			#add name to grpnms hash
			if ( exists( $grpnms{ $genotypex } ) ) {
				$grpnms{ $genotypex } = $grpnms{ $genotypex }."\t".$name;
			} else {
				$grpnms{ $genotypex } = $name;
			}

			#first deal with the genotype without taking sex into account
			if ( exists( $indivthrs{ $genotype } ) ) { #if there's already some data for that genotype
				$totals{ $genotype } = $totals{ $genotype } + 1; #add the new data
				$avgage{ $genotype } = $avgage{ $genotype } + $age{ $name }; #add on the age of this person
				my $ben;
				for $ben ( 0 .. $#thrsh ) {
					$indivthrs{ $genotype }[ $ben ] = $indivthrs{ $genotype }[ $ben ]."\t". $thrsh[ $ben ];
				}

			} else {                                      #else make the record
				$totals{ $genotype } = 1; 
				$avgage{ $genotype } = $age{ $name }; #get the age for this person
				$indivthrs{ $genotype } = [ @thrsh ]; #put the thresholds in as a new array ref
				$averages{ $genotype } = [ @thrsh ]; #get the averages array ready
				$stdevs{ $genotype } = [ @thrsh ]; #also get the standard deviation aray ready
			}

			#now repeat the process for the grouping by sex
			if ( exists( $indivthrs{ $genotypex } ) ) { #if there's already some data for that genotype
				$totals{ $genotypex } = $totals{ $genotypex } + 1; #add the new data
				$avgage{ $genotypex } = $avgage{ $genotypex } + $age{ $name }; #add on the age of this person
				my $ben;
				for $ben ( 0 .. $#thrsh ) {
					$indivthrs{ $genotypex }[ $ben ] = $indivthrs{ $genotypex }[ $ben ]."\t". $thrsh[ $ben ];
				}

			} else {                                      #else make the record
				$totals{ $genotypex } = 1; 
				$avgage{ $genotypex } = $age{ $name }; #get the age for this person
				$indivthrs{ $genotypex } = [ @thrsh ]; #put the thresholds in as a new array ref
				$averages{ $genotypex } = [ @thrsh ]; #get the averages array ready
				$stdevs{ $genotypex } = [ @thrsh ]; #also get the standard deviation array ready
			}
		}

		#now go through the totals and indivthrs hash
		foreach my $group ( keys %indivthrs ) {
			my $bunty;
			for $bunty ( 0 .. $#{ $indivthrs{ $group } } ) { #divide the sum of the thresholds by the total of that group, also divide by 10 and round to 2dp
				my $totalthresholds = 0; #count of valid thresholds, because some data have NAs
				my @thresholds = split(/\t/,$indivthrs{ $group }[ $bunty ]);
                                my $sumthr=0;
                                my $pip;
                                for $pip (0 .. $#thresholds) {
					if ( $thresholds[ $pip ] ne "NA" ) {
	                                        $thresholds[ $pip ] = $thresholds[ $pip ] /10; 
                                        	$sumthr = $sumthr+$thresholds[ $pip ];
						$totalthresholds = $totalthresholds + 1;
					}
                                }

				my $avg;
				if ( $indivthrs{ $group }[ $bunty ] !~ /\d/ ) { #catch the occasions where there's no digits, only NAs
					$avg = "NA";
					$totalthresholds = 1;
				} else {
	                               $avg = $sumthr/$totalthresholds; #get average
				}

				my $stdev = 0;
				my $sqtotal = 0;
				if ( $totalthresholds > 1 ) { #if there's more than one entry (otherwise standard deviation is 0)
					my $pin;
					for $pin ( 0 .. $#thresholds ) {
						if ( $thresholds[ $pin ] ne "NA" ) {
							$sqtotal = $sqtotal + (( $avg - $thresholds[ $pin ]) ** 2);
						}
					}
					$stdev = ($sqtotal/$totalthresholds) ** 0.5; 
				}

				if ( $avg ne "NA" ) {
					$avg = sprintf "%.2f", $avg;
				}
				$stdev = sprintf "%.2f", $stdev;
				$averages{ $group }[ $bunty ] = $avg;
				$stdevs{ $group }[ $bunty ] = $stdev;
				$indivthrs{ $group }[ $bunty ] = join("\t",@thresholds);

			}
			my $meanage = $avgage{ $group }/ $totals{ $group };
			$meanage = $meanage/100; #back to floating point
			$meanage = sprintf "%.2f", $meanage; #2dp
			$avgage{ $group } = $meanage; #store the average age
		}
	
		##deciding whether this is a candidate for outputting or not	
		my $printthis = 0;
		if ( ( $lookingforgene == 1 ) && ( $gene eq $genename )) { #if looking for a gene and this is the set gene, output everything
			$printthis= 1;
	               	foreach my $group ( keys %averages ) {
				if ( $group =~ /1/ ) {
					$notables{ $group } = 1; #setting the notables for every non-reference group so the individuals are output later
				}
			}
		} elsif ( ( $lookingforgene == 0 ) && ( $generategraphs == 1 )) { #if not looking for a specific gene but instead all variants should be output
			$printthis= 1;
	               	foreach my $group ( keys %averages ) {
				if ( $group =~ /1/ ) {
					$notables{ $group } = 1; #setting the notables for every non-reference group so the individuals are output later
				}
			}
		} elsif  ( ( $lookingforgene == 0 ) && ( exists( $averages{ "0_0" }) ) ) { #otherwise, cycle through again and do comparison, as long as there's a "wildtype" group (homozygous reference allele)
	               	foreach my $group ( keys %averages ) {
				if (( $group !~ /_N/ ) && ( $group =~ /1/ )) { #Only comparisons between the "wildtype" and an alternative allele group are wanted, in all participants together and groups with known sex
					my $counterright = 0; #counter for number of points passing criteria (right ear)
					my $counterleft = 0; #counter for number of points passing criteria (left ear)
					my @countthis = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0); #array to hold the indices of those points which pass the criteria for this group
       			                my $misty;
					for $misty ( 0 .. $#{ $averages{ $group } } ) {
						#by default, use the combined reference average
						my $avggrp = "0_0";
						#Use sex-specific reference averages instead if relevant and available - male and female only included at this time
						my $xchar = substr($group, -1); #get the final character of the group
						my $xgrp = "0_0_".$xchar; #create the relevant reference group based on that final character
						if ( exists( $averages{ $xgrp } )) { #if that group exists, use it
							$avggrp = $xgrp;
						}

						if ( ( $averages{ $group }[ $misty ] eq "NA" ) || ( $averages{ $avggrp }[ $misty ] eq "NA" ) ) { next; } #skip the stimuli for which there is only NA

						#if the difference in the thresholds for this stimulus is > set threshold and standard deviation of this group is < set sd, and there are at least 5 people in the group
						if ((( abs($averages{ $group }[ $misty ] - $averages{ $avggrp }[ $misty ] )) > $set ) && ($stdevs{ $group }[ $misty ] < $setstdv[ $misty ] ) && ( $totals{ $group } > 4 ) && ( $totals{ $avggrp } > 4 ) ) {

							if ( $averages{ $group} [ $misty ] > $averages{ $avggrp }[ $misty ] ) { #mark this as something to pay attention to - +1 for better hearing, -1 for worse
								$countthis[$misty ] = 1;
							} else {
								$countthis[$misty] = -1;
							}
							if ( $misty < 8 ) { #keep count of which ear these are happening in
								$counterright = $counterright + 1;
							} else {
								$counterleft = $counterleft + 1;
							}
						}
					}
					#if there are 2+ frequencies with different thresholds in each ear, add the current array of valid indices to the notables hash under their genotype
					if ( ( $counterright > 1 ) && ( $counterleft > 1 ) ) {
						$notables{ $group } = [ @countthis ];
						$printthis = 1; #variant passes filter and will be output
					}
				}
			}
		}

		if ( ( $printthis == 1 ) && ( $generategraphs == 0 ) ) { #variants will be printed in vcf format unless graphs are being generated
			print $outline;
		}
		if ( $generategraphs == 1 ) {
			my $filename = $data[0]."_".$data[1]."_".$data[2]."_".$data[3]."_".$data[4].".png";
			#output to text file
			print $varinfo."\t".$filename."\t".$genename."\n";
			foreach my $grpnym ( sort keys %grpnms ) {
				print "\t\t\t\t\t".$grpnym."\t".$grpnms{ $grpnym }."\n";
			}

			##Graphing
			my $multiChart = Chart::Gnuplot->new(
				output => "tmp.png",
				title => {
					text=>"$genename $data[2] $data[3]\>$data[4]",
					offset=>"0,-2",
					font=>"Helvetica, 40",
				},
				imagesize => '3,4',
			);
			my @charts=();
			#top left; left ear, groups
			$charts[0][0] =Chart::Gnuplot->new( 
				ylabel => {
					text => "Left ear threshold (dB HL)",
					offset => "-7,0",
					font => "Helvetica, 40",
				},
				legend => { 
					position => "bottom left",
					align => 'left',
					width => -6,
					height => 1,
					sample => {
						position => 'left',
					},
				},
				tmargin => 6,
				bmargin => 9,
			);
			#top right; right ear, groups
			$charts[0][1] =Chart::Gnuplot->new( 
				ylabel => {
					text => "Right ear threshold (dB HL)",
					font => "Helvetica, 40",
					offset => "-7,0",
				},
			);
			#bottom left - left ear individual thresholds
			$charts[1][0] =Chart::Gnuplot->new( 
				ylabel => {
					text => "Left ear threshold (dB HL)",
					offset => "-7,0",
					font => "Helvetica, 40",
				},
			);
			#bottom right; right ear, individual thresholds
			$charts[1][1] =Chart::Gnuplot->new( 
				ylabel => {
					text => "Right ear threshold (dB HL)",
					font => "Helvetica, 40",
					offset => "-7,0",
				},
			);

			#format the charts
			for my $xchart ( 0 .. 1 ) {
				for my $ychart ( 0 .. 1 ) {
					$charts[ $xchart ][ $ychart ]->set(
						xlabel => {
							text => "Frequency (kHz)",
							font => "Helvetica, 40",
							offset => "0, -3",
						},
						grid => {
							linetype => 1,
							width => 2,
							color => '#404040',
						},
						xrange => [0.5, $xmx],
						yrange => [135, -20],
						ytics => {
							labels => ['"120" 120','"" 110','"100" 100','"" 90','"80" 80','"" 70','"60" 60','"" 50','"40" 40','"" 30','"20" 20','"" 10','"0" 0','"" -10'],
							fontcolor => "black",
							length => "0",
							font => "Helvetica, 40",
						},
						bars => 4, 
						size =>'square 0.5',
						border => {
							color => '#404040',
							width => 3,
						},
						tmargin => 6,
						bmargin => 9,
					);
					$charts[$xchart][$ychart]->line( from => "0.5,0", to=> $maxx, width=>3, color=>'#000000',);
				}
			}

                        #set the x axis labels
                        for my $xchaxis ( 0 .. 1 ) {
                                for my $ychaxis ( 0 .. 1 ) {
                                        if ( $start125 ) {
                                                $charts[$xchaxis][$ychaxis]->set(
                                                        xtics => {
                                                                labels => ['".125" 1','".25" 2','".5" 3','"1" 4','"2" 5','"4" 6','"6" 6.56','"8" 7'],
                                                                offset => "0,-1",
                                                                fontcolor => "black",
                                                                font => "Helvetica, 40",
                                                                length => "0",
                                                        },
                                                );
                                        } else {
                                                $charts[$xchaxis][$ychaxis]->set(
                                                        xtics => {
                                                                labels => ['".25" 1','".5" 2','"1" 3','"2" 4','"3" 4.56','"4" 5','"6" 5.56','"8" 6'],
                                                                offset => "0,-1",
                                                                fontcolor => "black",
                                                                font => "Helvetica, 40",
                                                                length => "0",
                                                        },
                                                );
                                        }
                                }
                        }


			#go through ordering array and pull out the existing groups in the desired order
			my @order;
			my $pod;
			for $pod ( 0 .. $#ordering ) {
				if ( exists( $averages{ $ordering[ $pod ] } ) ) { push(@order, $ordering[ $pod ]); }
			}
			
			my $ycoord = -10; #y coordinate for adding significance stars, goes up by 2 every time a group has stars added so they all should be visible

			my $pip;
			for $pip ( 0 .. $#order) {
				my $group = $order[ $pip ]; 
				my @avgs = @{$averages{ $group }};
				my @stdvs = @{$stdevs{ $group }};
				my @thrsh = @{$indivthrs{ $group }};
				my $n = $totals{ $group };

				if ( ( $group =~ /1/ ) || ( $group =~ /0_0/ ) ) { #if the group is the reference or has an alt allele
					my $colour="#F7F056"; #default
					my $pointtype="square"; #default
					my $gtype=$group;
					if ( exists( $colours{ $group } ) ) {
						$colour=$colours{$group};
					}
					if ( exists( $gtypes{ $group } ) ) {
						if ( $ismit ) {
							$gtype=$gtypes{"MT".$group};
						} else {
							$gtype=$gtypes{$group};
						}
					}
					if ( exists( $points{ $group } ) ) {
						$pointtype=$points{$group};
					}

					#adjust the average age to use a range if n=1, so that individual ages are not displayed
                                        
					my $legage = $avgage{ $group };
					if ( $n == 1 ) {
						if ( $legage < 40 ) { $legage = "<40"; }
						elsif ( $legage < 45 ) { $legage = "40-45"; }
						elsif ( $legage < 50 ) { $legage = "45-50"; }
						elsif ( $legage < 55 ) { $legage = "50-55"; }
						elsif ( $legage < 60 ) { $legage = "55-60"; }
						elsif ( $legage < 65 ) { $legage = "60-65"; }
						elsif ( $legage < 70 ) { $legage = "65-70"; }
						elsif ( $legage < 75 ) { $legage = "70-75"; }
						elsif ( $legage < 80 ) { $legage = "75-80"; }
						elsif ( $legage < 85 ) { $legage = "80-85"; }
						elsif ( $legage < 90 ) { $legage = "85-90"; }
						else { $legage= ">90"; }
					}

					my $legtitle = $gtype." n=".$n.", ".$legage."y";

					my $sig="*";

					if ( $group =~ /F/ ) {
						$sig="=";
					} elsif ( $group =~ /M/ ) {
						$sig = "+";
					}

					#if this is one of the groups with "notable" thresholds, add stars to top of graph to denote which ones count
					#don't do this if "generategraphs" is set
					if ( (exists( $notables{ $group }) ) && ( $generategraphs == 0) ) {
						my $sab;
						for $sab ( 0 .. 7 ) {
							if ( $notables{ $group }[ $sab ] != 0 ) {
								$charts[0][1]->label(
									text => $sig,
									fontcolor => $colour,
									position=> "$x[$sab], $ycoord center", 
								)
							}
						}
						my $bas;
						for $bas ( 8 .. $#{ $notables{ $group } } ) {
							if ( $notables{ $group }[ $bas ] != 0 ) {
								my $cind = $bas-8;
								$charts[0][0]->label(
									text => $sig,
									fontcolor => $colour,
									position=> "$x[$cind], $ycoord center", 
								)
							}
						}
						$ycoord = $ycoord - 2;
					}

					#Although Gnuplot handles internal missing NA values fine, it can't cope with them at the start of a dataset so these have to be checked
					my @g_rdata;
					my @g_ldata;
					my @gr_indx;
					my @gl_indx;
					my @grerr;
					my @glerr;

					#go through data and catch starting NAs
					for my $grpnts ( 0 .. 7) {
						if ( $avgs[ $grpnts ] ne "NA" ) {
							push( @g_rdata, $avgs[ $grpnts ]);
							push( @grerr, $stdvs[ $grpnts ]);
							push( @gr_indx, $x[ $grpnts ]);
						}
					}
					for my $glpnts ( 8 .. $#avgs ) {
						if ( $avgs[ $glpnts ] ne "NA" ) {
							push( @g_ldata, $avgs[ $glpnts ]);
							push( @glerr, $stdvs[ $glpnts ]);
							push( @gl_indx, $x[ ($glpnts - 8) ]);
						}
					}

					my @left;
					my @right;
					@left[0] = \@g_ldata;
					@right[0] = \@g_rdata;
					@left[1] = \@glerr;
					@right[1] = \@grerr;			

					my $dataset_l = Chart::Gnuplot::DataSet->new (
						xdata => \@gl_indx,
						ydata => \@left,
						style => "yerrorlines",
						width => 7,
						color => $colour,
						pointtype => $pointtype,
						pointsize => 4,
						title => $legtitle,
					);
					$charts[0][0] -> add2d($dataset_l);

					my $dataset_r = Chart::Gnuplot::DataSet->new (
						xdata => \@gr_indx,
						ydata => \@right,
						style => "yerrorlines",
						width => 7,
						color => $colour,
						pointtype => $pointtype,
						pointsize => 4,
					);
					$charts[0][1] -> add2d($dataset_r); 

					if ( $group =~ /0_0/ ) { #if it's a reference group, plot on individual chart as well (no legend applied this time)
						my $dataset_2l = Chart::Gnuplot::DataSet->new (
	                                                xdata => \@gl_indx,
        	                                        ydata => \@left,
                	                                style => "yerrorlines",
                        	                        width => 7,
                                	                color => $colour,
                                        	        pointtype => $pointtype,
                                                	pointsize => 4,
	                                        );

						$charts[1][0] -> add2d($dataset_2l);
						$charts[1][1] -> add2d($dataset_r);
					}
				}
			}

			#now go through the individuals carrying non-reference alleles and plot
			foreach my $grpkey ( keys %notables ) { #for each group
				my @peeps = ();
				my $combined=0;
				if ( $grpkey !~ /M|F/ ) { #if this is one of the combined sex groups, get the individual names from both
					my $keym = $grpkey."_M";
					my $keyf = $grpkey."_F";
					my $keyn = $grpkey."_N";
					if (exists( $grpnms{ $keyf } )) {
						my ( @pps2 ) = split(/\t/, $grpnms{ $keyf } );
						push( @peeps, @pps2 );
					}
					if (exists( $grpnms{ $keym } )) {
						my ( @pps2 ) = split(/\t/, $grpnms{ $keym } );
						push( @peeps, @pps2 );
					}
					if (exists( $grpnms{ $keyn } )) {
						my ( @pps2 ) = split(/\t/, $grpnms{ $keyn } );
						push( @peeps, @pps2 );
					}
					$combined = 1; #flag that individuals need to be assigned point types by sex
				} else {
					(@peeps) = split(/\t/, $grpnms{ $grpkey }); #get the people in that group
				}
				for my $faff ( 0 .. $#peeps ) { #go through the people in that group
					my @pdata = @{$originalthrsh{ $peeps[ $faff ] }}; #get their original data
					my $indcol = $colours{ $grpkey };
					my $indpt = $points{ $grpkey }; 
					if ( $combined ) {
						my $indkey = $grpkey."_".$sex{ $peeps[ $faff] };
						$indcol = $colours{ $indkey };
						$indpt = $points{ $indkey };
					}

					#Although Gnuplot handles internal missing NA values fine, it can't cope with them at the start of a dataset so these have to be checked
					my @p_rdata;
					my @p_ldata;
					my @r_indx;
					my @l_indx;

					#go through data and catch starting NAs
					for my $rpnts ( 0 .. 7) {
						if ( $pdata[ $rpnts ] ne "NA" ) {
							push( @p_rdata, $pdata[ $rpnts ]);
							push( @r_indx, $x[ $rpnts ]);
						}
					}
					for my $lpnts ( 8 .. $#pdata ) {
						if ( $pdata[ $lpnts ] ne "NA" ) {
							push( @p_ldata, $pdata[ $lpnts ]);
							push( @l_indx, $x[ ($lpnts - 8) ]);
						}
					}
					
					my $dataset_indl = Chart::Gnuplot::DataSet->new (
						xdata => \@l_indx,
						ydata => \@p_ldata,
						style => "linespoints",
						width => 7,
						color => $indcol,
						pointtype => $indpt,
						pointsize => 4,
					);
					my $dataset_indr = Chart::Gnuplot::DataSet->new (
						xdata => \@r_indx,
						ydata => \@p_rdata,
						style => "linespoints",
						width => 7,
						color => $indcol,
						pointtype => $indpt,
						pointsize => 4,
					);
					$charts[1][0] -> add2d($dataset_indl);
					$charts[1][1] -> add2d($dataset_indr);
				}
			}

			$multiChart->multiplot(\@charts);

			#imagemagick call to flatten output image
			my $cmd = "convert tmp.png -background white -flatten ".$filename;
			system($cmd);
		}
	}
}
$fhc->close();
