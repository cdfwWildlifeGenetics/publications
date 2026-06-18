# NCBI submission reordering scaffolds

**NOTES: splitting scaffolds to have FCS remove and trim contamination**

**Description: write each scaffold out to its own file.  This will allow me to edit the scaffolds I need to.  I will then concat all the scaffold files together for NCBI submission.**

*Load genome assembly to NCBI. Instructions: https://submit.ncbi.nlm.nih.gov/subs/genome/ I used web browser option. I submitted the genome and waited for NCBI's pipeline FCS to print out an error to determine which scaffolds I had to edit.  If there is contamination in the center of a scaffold then you will have to split that scaffold, create a new fasta file of the genome, and resubmit.  If contamination is at the ends of the scaffold FCS will trim these regions.*


## Step 1:

```
Trim:
Sequence name, length, span(s), apparent source
PGA_scaffold_16__97_contigs__length_84625481	84625481	83810465..83882985	plnt:plants
PGA_scaffold_17__81_contigs__length_66046125	66046125	63929032..63929074	adaptor:NGB01029.1-not_cleaned
PGA_scaffold_17__81_contigs__length_66046125	66046125	65524839..65709940	plnt:plants
PGA_scaffold_20__217_contigs__length_76918119	76918119	75540516..75564286	plnt:plants
PGA_scaffold_22__119_contigs__length_66889936	66889936	34354253..34360176	plnt:plants
PGA_scaffold_29__303_contigs__length_46869212	46869212	45506421..45513568	plnt:plants
PGA_scaffold_8__74_contigs__length_47778798	47778798	17774712..17774739	adaptor:NGB04004.1-not_cleaned
```
___

## Step 2: install / load packages
```
source /share/cdfwwildlife/hallas_dedicated/Miniconda/etc/profile.d/conda.sh

conda activate popgen

module load seqtk
```
___

## Step 3: Create a list of scaffold names using bioawk
```
bioawk -c fastx '{ print $name }' /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa > /share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/scaffold_list.txt
```
___


## Step 4: Create individual scaffold files
```
while read scaffold; do
    seqtk subseq /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa <(echo "$scaffold") \
    > "/share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/${scaffold}.fa"
done < /share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/scaffold_list.txt
```
___

## Step 5: Split scaffolds at contamination region
**NOTES: This could probably be done more efficient but I ran a separate command for each scaffold that needed to be split.** 

*PGA_scaffold_16__97_contigs__length_84625481	84625481	83810465..83882985*
```
bioawk -c fastx '$name == "PGA_scaffold_16__97_contigs__length_84625481" { print ">"$name"_part1\n"substr($seq, 1, 83810464) }' /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa > PGA_scaffold_16_part1.fa
bioawk -c fastx '$name == "PGA_scaffold_16__97_contigs__length_84625481" { print ">"$name"_part2\n"substr($seq, 83810465, 84625481) }' /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa > PGA_scaffold_16_part2.fa
```
___

*PGA_scaffold_17__81_contigs__length_66046125	66046125	63929032..63929074*
```
bioawk -c fastx '$name == "PGA_scaffold_17__81_contigs__length_66046125" { print ">"$name"_part1\n"substr($seq, 1, 63929031) }' /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa > PGA_scaffold_17_part1.fa
bioawk -c fastx '$name == "PGA_scaffold_17__81_contigs__length_66046125" { print ">"$name"_part2\n"substr($seq, 63929032, 66046125) }' /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa > PGA_scaffold_17_part2.fa
```
___

*PGA_scaffold_20__217_contigs__length_76918119	76918119	75540516..75564286*
```
bioawk -c fastx '$name == "PGA_scaffold_20__217_contigs__length_76918119" { print ">"$name"_part1\n"substr($seq, 1, 75540515) }' /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa > PGA_scaffold_20_part1.fa
bioawk -c fastx '$name == "PGA_scaffold_20__217_contigs__length_76918119" { print ">"$name"_part2\n"substr($seq, 75540516, 76918119) }' /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa > PGA_scaffold_20_part2.fa
```
___

*PGA_scaffold_22__119_contigs__length_66889936	66889936	34354253..34360176*
```
bioawk -c fastx '$name == "PGA_scaffold_22__119_contigs__length_66889936" { print ">"$name"_part1\n"substr($seq, 1, 34354252) }' /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa > PGA_scaffold_22_part1.fa
bioawk -c fastx '$name == "PGA_scaffold_22__119_contigs__length_66889936" { print ">"$name"_part2\n"substr($seq, 34354253, 66889936) }' /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa > PGA_scaffold_22_part2.fa
```
___

*PGA_scaffold_29__303_contigs__length_46869212	46869212	45506421..45513568*
```
bioawk -c fastx '$name == "PGA_scaffold_29__303_contigs__length_46869212" { print ">"$name"_part1\n"substr($seq, 1, 45506420) }' /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa > PGA_scaffold_29_part1.fa
bioawk -c fastx '$name == "PGA_scaffold_29__303_contigs__length_46869212" { print ">"$name"_part2\n"substr($seq, 45506421, 46869212) }' /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa > PGA_scaffold_29_part2.fa
```
___

*PGA_scaffold_8__74_contigs__length_47778798	47778798	17774712..17774739*
```
bioawk -c fastx '$name == "PGA_scaffold_8__74_contigs__length_47778798" { print ">"$name"_part1\n"substr($seq, 1, 17774711) }' /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa > PGA_scaffold_8_part1.fa
bioawk -c fastx '$name == "PGA_scaffold_8__74_contigs__length_47778798" { print ">"$name"_part2\n"substr($seq, 17774712, 47778798) }' /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa > PGA_scaffold_8_part2.fa
```
___

## Step 6: Remove unsplit scaffolds
```
rm PGA_scaffold_16__97_contigs__length_84625481.fa
rm PGA_scaffold_17__81_contigs__length_66046125.fa
rm PGA_scaffold_20__217_contigs__length_76918119.fa
rm PGA_scaffold_22__119_contigs__length_66889936.fa
rm PGA_scaffold_29__303_contigs__length_46869212.fa
rm PGA_scaffold_8__74_contigs__length_47778798.fa
```
___

## Step 7: Creates list of File, Scaffold name, and length
```
for file in /share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/*.fa; do
    bioawk -c fastx -v fname=$(basename "$file") '{ print fname "\t" $name "\t" length($seq) }' "$file" >> scaffold_lengths.tsv
done
```
___

## Step 8: Puts scaffolds in descending order*
```
sort -k3,3nr /share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/scaffold_lengths.tsv -o /share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/scaffold_lengths.tsv

```
___

## Step 9: Concat scaffolds
*Concats scaffolds based on scaffold length file.  This keeps scaffolds in descending order*
```
awk '{print $1}' /share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/scaffold_lengths.tsv | while read -r filename; do
    cat "/share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/$filename" >> /share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/CDFW_2019_HSUGinger_reordered.fa
done
```
___

## Step 10: Renames scaffolds based on order and add length*
```
bioawk -c fastx '{ print ">scaffold-" ++i"_"length($seq)"\n"$seq }' < /share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/CDFW_2019_HSUGinger_reordered.fa > /share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/CDFW_2019_HSUGinger_new.fa
```
___

## Step 11: check lengths. they should be equal*
```
bioawk -c fastx '{print length($seq)}' /share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/CDFW_2019_HSUGinger_reordered.fa | awk '{sum += $1} END {print sum}'

bioawk -c fastx '{print length($seq)}' /share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/CDFW_2019_HSUGinger_new.fa | awk '{sum += $1} END {print sum}'

bioawk -c fastx '{print length($seq)}' /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa | awk '{sum += $1} END {print sum}'
```
___

## Step 12: remove scaffold files
```
rm PGA_scaffold_*
rm CDFW_2019_HSUGinger_reordered.fa
```
___

## Step 13: move genome to mac
```
scp barbera:/share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/CDFW_2019_HSUGinger.fa /Users/jhallas/Documents/cdfw/CCGP_genome/NCBI_submission 
```
___


# third Submission

**IMPORTANT NOTES: after resubmission to NCBI there was still a scaffold that had contamination in the middle of it.  Following script fixes this issue**

```
Trim:
Sequence name, length, span(s), apparent source
scaffold-43	2117051	1595765..1780866	plnt:plants
```

**create file with only scaffold-43.  this file will be split at contamination interval**
```
echo "scaffold-43" > scaffold-43.list

seqtk subseq CDFW_2019_HSUGinger.fa scaffold-43.list > scaffold-43.fa
```

**create a list of scaffolds excluding scaffold-43**
```
bioawk -c fastx '$name != "scaffold-43" { print $name }' CDFW_2019_HSUGinger.fa > CDFW_2019_HSUGinger_V2keep.list
```

**subsets genome with only scaffold on the keep list**
```
seqtk subseq CDFW_2019_HSUGinger.fa CDFW_2019_HSUGinger_V2keep.list > CDFW_2019_HSUGingerV2.fa
```

**split scaffold-43.fa into two separate files**
```
bioawk -c fastx '$name == "scaffold-43" { print ">"$name"_part1\n"substr($seq, 1, 1595764) }' CDFW_2019_HSUGinger.fa > scaffold-43_part1.fa
bioawk -c fastx '$name == "scaffold-43" { print ">"$name"_part2\n"substr($seq, 1595765, 2117051) }' CDFW_2019_HSUGinger.fa > scaffold-43_part2.fa

```

**concat scaffold-43 part1 and part2 to CDFW_2019_HSUGingerV2.fa**
```
cat scaffold-43_part1.fa >> CDFW_2019_HSUGingerV2.fa
cat scaffold-43_part2.fa >> CDFW_2019_HSUGingerV2.fa
```

**Reorder and label scaffolds in CDFW_2019_HSUGingerV2.fa**
```
bioawk -c fastx '{ print length($seq), $name }' CDFW_2019_HSUGingerV2.fa | sort -nr | awk '{print $2}' > CDFW_2019_HSUGingerV2_concat_sorted_names.list

seqtk subseq CDFW_2019_HSUGingerV2.fa CDFW_2019_HSUGingerV2_concat_sorted_names.list > CDFW_2019_HSUGingerV2_sorted_output.fa
```

**rename scaffolds**
```
bioawk -c fastx '{ print ">scaffold-" ++i"_"length($seq)"\n"$seq }' < /share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/CDFW_2019_HSUGingerV2_sorted_output.fa > /share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/CDFW_2019_HSUGingerV3.fa
```

**check lengths**
```
bioawk -c fastx '{print length($seq)}' /share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/CDFW_2019_HSUGingerV3.fa | awk '{sum += $1} END {print sum}'

bioawk -c fastx '{print length($seq)}' /share/cdfwwildlife/Deer_Genome_Assembly/reference.filtered.fa | awk '{sum += $1} END {print sum}'
```

**transfer to mac**
```
scp barbera:/share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/CDFW_2019_HSUGingerV3.fa /Users/jhallas/Documents/cdfw/CCGP_genome/NCBI_submission 
```

*clean up directory after submission to NCBI!!!!*
```
rm scaffold-43.fa
rm scaffold-43_part1.fa
rm scaffold-43_part2.fa
rm CDFW_2019_HSUGingerV2_sorted_output.fa
```
___

**Resubmit new genome**
___

# Upload to NCBI SRA FASTQ files
**NOTE: this was annoying.  also took a lot of time so make sure your computer doesn't turn off and end the upload. follow instructions on https://submit.ncbi.nlm.nih.gov/subs/sra/SUB15279929/files**

*1-download "Download now for macOS".  I downloaded this to my personal account.  The code used for NCBI submission is not in the path.  need to call command using full path.  It uses spaces which is stupid.*


*download fastq to personal computer*
```
scp barbera:/share/cdfwwildlife/Deer_Genome_Assembly/00_raw_sequencing_data/10X_reads/Black-tailed_deer_S2_L005_R1_001.fastq.gz /Users/jhallas/Documents/cdfw/CCGP_genome/NCBI_submission/Fastq 


scp barbera:/share/cdfwwildlife/Deer_Genome_Assembly/00_raw_sequencing_data/10X_reads/Black-tailed_deer_S2_L005_R2_001.fastq.gz /Users/jhallas/Documents/cdfw/CCGP_genome/NCBI_submission/Fastq 
```

*step 1: download Aspera Connect Software*
*step 2: download key file*
*step 3: You may use the following command to upload files via Aspera Command-Line*
```
/Users/jhallas/Applications/IBM\ Aspera\ Connect.app/Contents/Resources/ascp -i /Users/jhallas/Documents/cdfw/CCGP_genome/NCBI_submission/key/aspera.openssh -QT -l100m -k1 -d /Users/jhallas/Documents/cdfw/CCGP_genome/NCBI_submission/fastq2 subasp@upload.ncbi.nlm.nih.gov:uploads/joshua.hallas_gmail.com_77XLYxhq
```
___

# QUAST to evaluate genomes used by Merly
*Install*
```
source /share/cdfwwildlife/hallas_dedicated/Miniconda/etc/profile.d/conda.sh

conda activate popgen

conda install bioconda::quast
conda install bioconda::purge_dups
```
___
*count scaffold lengths and remove mtDNA scaffold*
```
module load seqtk

bioawk -c fastx '{ print $name, length($seq) }' /share/cdfwwildlife/Deer_CCGP_Genome_Assembly/GCA_042768445.1_mOdoHem1.0.p_genomic.fna  | sort --key 2 --numeric-sort --reverse | awk '{sum+=$1; lengths[NR]=$1} END {half=sum/2; for (i=1; i<=NR; i++) {cumsum+=lengths[i]; if (cumsum >= half) {print lengths[i]; break}}}'

bioawk -c fastx '{ print $name }' GCA_042768445.1_mOdoHem1.0.p_genomic.fna | sed 's/CM091143.1//g' > scaffold_names.list

seqtk subseq GCA_042768445.1_mOdoHem1.0.p_genomic.fna scaffold_names.list > GCA_042768445.1_mOdoHem1.0.p_genomicNOMito.fna
```

*Run quast with estimated genome size.  Used estimated ref length based on value used by Merly*
```
quast --est-ref-size 2980000000 -o /share/cdfwwildlife/hallas_dedicated/ccgp_genome/quast_output/quast_output_all_5June2025 /share/cdfwwildlife/Deer_CCGP_Genome_Assembly/GCA_042768445.1_mOdoHem1.0.p_genomicNOMito.fna /share/cdfwwildlife/Deer_Genome_Assembly/NCBI_submission/CDFW_2019_HSUGinger00000000.fsa /share/cdfwwildlife/hallas_dedicated/reference_genomes/Lamb_et_al_2021_ref/finalmuledeergenome_filtered.fasta /share/cdfwwildlife/hallas_dedicated/reference_genomes/Russel_et_al_2019/GCA_004115125.1_UofA_Ohem_1.0_genomic.fna /share/cdfwwildlife/hallas_dedicated/reference_genomes/Sitka_PRJNA476345/GCA_003697985.1_ASM369798v1_genomic.fna /share/cdfwwildlife/hallas_dedicated/reference_genomes/London_et_al_2022/GCF_023699985.2_Ovbor_1.2_genomic.fna /share/cdfwwildlife/hallas_dedicated/reference_genomes/Seabury_et_al_2011_ref/GCF_002102435.1_Ovir.te_1.0_genomic.fna

scp -r jhallas@barbera.genomecenter.ucdavis.edu:/share/cdfwwildlife/hallas_dedicated/ccgp_genome/quast_output/quast_output_all_5June2025 /drives/o/Genetics\ Lab/Wildlife\ Genetics\ Program/Josh/Deer_CCGP_genome
```
___
