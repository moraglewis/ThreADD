# ThreADD

These scripts were written to compare the average pure tone thresholds of people separated into groups by genotype and sex for a list of input variants. Required input files include the pure tone thresholds and the variant calls. 

The first script (ThreADD.pl) does the comparisons and outputs a list of variants which pass the user-supplied minimum threshold and maximum standard deviation (adjusted by frequency). It can also output graphs without doing comparisons, in which case it will either output a graph for every input variant, or it will output a graph for every input variant within a specific gene.

The second script (ThreADD_permute.pl) carries out the same comparisons followed by a permutation test. Only variants which pass both the comparison test and the permutation test are output as graphs.

![threadd](https://user-images.githubusercontent.com/29064421/216606662-2bc22822-ea6c-4e36-a7a5-ff4bad661e22.png)

The diagram is taken from a poster to be presented at the Association for Research in Otolaryngology Midwinter Meeting 2023. A paper describing the work is in progress.

ThreADD.pl takes as input:
1. A tab-separated file containing thresholds in the following format:
   - A header line starting "SEX" and containing sex, person ID, secondary ID, collection date, age at collection, noise history, right ear thresholds (for 8 pure tones in Hz), space, left ear thresholds (for the same 8 pure tones). 
   - Values for the SEX column can be anything, but only "M" and "F" values will be used when assigning individuals to sex-specific groups. All IDs which are neither "M" nor "F" will be assigned to the "N" group, which is not used for sex-specific comparison.
   - Any unavailable threshold values should be "NA". Any other unavailable values should be ".". The script only uses sex, age and ID in addition to the tone thresholds.
   - The default 8 tones are 0.25, 0.5, 1, 2, 3, 4, 6 and 8kHz. If the first tone is 125Hz, then the script switches to use 0.125, 0.25, 0.5, 1, 2, 4, 6 and 8kHz. No other options are provided.
2. A vcf file with variants. The IDs in the vcf header should correspond to the IDs in the threshold file, but it is not a problem if there are extra IDs in either file. 
   - Note that this script does not handle multiallelic sites, but vcf files can be "flattened" like so:
     - G   A,C     0/1     0/2     1/2
     - becomes:
     - G   A       0/1     0/v     1/v
     - G   C       0/v     0/1     v/1
     - "v" here means any allele which is neither the reference nor the alternate
   - "Minority calls" (referred to in comments) are calls which are not simple ones, eg (0/1, 1/1), and indicate a failed call
   -  The script will handle mitochondrial variants, but expects them to be marked "HET1", "HOM1", "WT" or "FAIL" (for heteroplasmic, homoplasmic, wildtype or fail)
   - If variants are plotted, a gene name from the VCF file will be used as a title for the plot. A VEP-annotated vcf file will work, or you can provide your own names in the INFO field. The script expects values to be separated by "|", and the gene name and Ensembl ID should be in the fourth and fifth entry.
     - For example:
     - 1       1485777 rs1622213       G       A       100.5 PASS    CSQ=A|splice_region_variant|LOW|ATAD3B|ENSG00000160072
3. A desired minimum threshold difference, in dB
4. A maximum permitted standard deviation. This is adjusted by the script to allow for more variability at high frequencies (+10) and less at low frequencies (-5).
5. (Optional) A gene name (eg ATAD3B) or "gengr". 
   - In the case of a gene name being input, the script will generate graphs for every variant in that gene in the file.
   - In the case of "gengr", the script will generate graphs for all variants in the file. No comparisons will be carried out. 
   - In both cases, the set threshold and standard deviation still need to be entered, but will not be used.

The output is a list of variants which pass the required settings, and can be used as input for ThreADD_permute.pl. If a gene name or "gengr" is supplied as the fifth argument, image files with the graphs in will also be output.
Please note that ImageMagick is required for the image generation.
If graphs are generated without testing thresholds (eg by inputting a gene name or "gengr" as the fifth argument), the output text file instead contains details of each variant for which a graph was generated, its filename and the participant IDs in each grouping ("0_1_F" means female participants with a 0/1 genotype, ie heterozygotes).
If any group is plotted which contains only one participant, the "average" age of that group is not reported to protect privacy. Instead, a range is reported (eg 67 years becomes 65-70y, and 80 years becomes 80-85y).

ThreADD_permute.pl takes as input:
1. A tab-separated file containing thresholds in the following format:
   - A header line starting "SEX" and containing sex, person ID, secondary ID, collection date, age at collection, noise history, right ear thresholds (for 8 pure tones in Hz), space, left ear thresholds (for the same 8 pure tones)
   - Values for the SEX column can be anything, but only "M" and "F" values will be used when assigning individuals to sex-specific groups. All IDs which are neither "M" nor "F" will be assigned to the "N" group, which is not used for sex-specific comparison.
   - Any unavailable threshold values should be "NA". Any other unavailable values should be ".". The script only uses sex, age and ID in addition to the tone thresholds.
   - The default 8 tones are 0.25, 0.5, 1, 2, 3, 4, 6 and 8kHz. If the first tone is 125Hz, then the script switches to use 0.125, 0.25, 0.5, 1, 2, 4, 6 and 8kHz. No other options are provided.
2. A vcf file with variants. The IDs in the vcf header should correspond to the IDs in the threshold file, but it is not a problem if there are extra IDs in either file. 
   - Note that this script does not handle multiallelic sites, but vcf files can be "flattened" like so:
     - G   A,C     0/1     0/2     1/2
     - becomes:
     - G   A       0/1     0/v     1/v
     - G   C       0/v     0/1     v/1
     - "v" here means any allele which is neither the reference nor the alternate
   - "Minority calls" (referred to in comments) are calls which are not simple ones, eg (0/1, 1/1), and indicate a failed call
   -  The script will handle mitochondrial variants, but expects them to be marked "HET1", "HOM1", "WT" or "FAIL" (for heteroplasmic, homoplasmic, wildtype or fail)
   - If variants are plotted, a gene name from the VCF file will be used as a title for the plot. A VEP-annotated vcf file will work, or you can provide your own names in the INFO field. The script expects values to be separated by "|", and the gene name and Ensembl ID should be in the fourth and fifth entry.
     - For example:
     - 1       1485777 rs1622213       G       A       100.5 PASS    CSQ=A|splice_region_variant|LOW|ATAD3B|ENSG00000160072
3. A desired minimum threshold difference, in dB
4. A maximum permitted standard deviation. This is adjusted by the script to allow for more variability at high frequencies (+10) and less at low frequencies (-5).
5. A file of random numbers, one per line, with no duplicates. This is used to seed the randomisation of the assignment of individuals to groups in the permutation, and the number of lines in the file dictates how many permutations are carried out.

The output is a text file containing details of the variants which passed the permutation filter, the filename of the graph generated, the number of permutations which resulted in a similar average threshold difference, the gene name and the participant IDs in each grouping ("0_1_F" means female participants with a 0/1 genotype, ie heterozygotes).
If any group is plotted which contains only one participant, the "average" age of that group is not reported to protect privacy. Instead, a range is reported (eg 67 years becomes 65-70y, and 80 years becomes 80-85y).
Also output are image files, one per filtered variant, showing the audiograms with averages, and also showing individual thresholds for those groups which were sufficiently different to pass the initial filter (the set threshold and standard deviation).

