# Run rehh on polarized and BEAGLE-phased data
Scripts:
- **polarize_alleles.slurm** - polarize alleles based on PRE major allele
- **runREHH.slurm** - masster script; submit beagle & rehh scripts with dependencies to slurm
- **runBEAGLE.sbatch** - phase genotypes with beagle
- **REHH_iHH.R** - calculate iHH per population per scaffold
- **REHH_iHH.sbatch** - execute **REHH_iHH.R** as an array of jobs
- **REHH_iHS.R** - combine iHH across scaffolds for each population, calculate iHS, and find per pop candidate regions
- **REHH_iHS.sbatch** - execute **REHH_iHS.R**
- **REHH_Rsb.R** - read in iHS for each population and calculate pariwise Rsb & XP-EHH
- **REHH_Rsb.sbatch** - execute **REHH_Rsb.R**
- **REHH_Rsb_filt.sbatch** - filter 10Kb smoothed Rsb & XP-CLR output to 10% overlapping windows

## Polarize alleles based on PRE major allele - **polarize_alleles.slurm**
```bash
#!/bin/bash

#SBATCH --job-name=polarize
#SBATCH --nodes=1
#SBATCH -t 10-00:00:00
#SBATCH --mem=50GB
#SBATCH --output=/share/cdfwwildlife/MYLU_NovaSeq/SCRIPTS/slurmout/polarize_%A.out
#SBATCH --error=/share/cdfwwildlife/MYLU_NovaSeq/SCRIPTS/slurmout/polarize_%A.err

start=`date +%s`
echo $HOSTNAME
dir="/share/cdfwwildlife/MYLU_NovaSeq"
cnd="/share/cdfwwildlife/Capel_Dedicated/miniconda"

aklog
source $cnd/etc/profile.d/conda.sh
conda activate osjdk17
module load bcftools
module load vcftools

## find major allele for PRE pops
bcftools view -S $dir/Sample_Lists/PRE.ind -O z gatk.snp.qual_hard_filtered_autosomes.vcf.gz | \
    vcftools --gzvcf - --freq --out $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE
tail -n +2 $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE.frq | \
    tr ":" "\t" | cut -f 1,2,5- | \
    awk '{if ($4 < 0.5) print $1"\t"$2"\t"$5"\t"$3; else print $1"\t"$2"\t"$3"\t"$5}' | \
    awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | bgzip > $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.bed.gz
tabix $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.bed.gz

## polarize REF allele in VCF
echo "Adding AA INFO field to vcf..."
bcftools view $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.vcf.gz | \
    vcf-annotate -a $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.bed.gz \
                 -d key=INFO,ID=AA,Number=1,Type=String,Description='Ancestral Allele' \
                 -c CHROM,FROM,TO,INFO/AA | \
    bcftools view -O z -o $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.vcf.gz

echo "Replacing REF and ALT alleles with AA..."
$cnd/envs/osjdk17/bin/java -jar $dir/SCRIPTS/JVARKIT/jvarkit.jar vcffilterjdk \
                           -f $dir/SCRIPTS/test_ancestral_allele/script.js $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.vcf.gz | \
    bgzip > $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.polarized.vcf.gz
rm $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.vcf.gz

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime"
```

## Phase genotypes by scaffold & population - **runBEAGLE.sbatch**
Some notes on running beagle v5.4:
- Increasing the number of iterations and the number of phase states increases accuracy, but also memory usage
- Reducing the window size reduces the memory usage - reduce to be able to increase # iterations & phase states
- Using more threads will decrease run time, but can also cause a Java memory usage error (heap space runs out) - limit # threads
```bash
#SBATCH --job-name=beagle
#SBATCH --mem=10G
#SBATCH --array=1-168
#SBATCH --ntasks-per-node=8

dir="/share/cdfwwildlife/MYLU_NovaSeq"
ar=`sed "${SLURM_ARRAY_TASK_ID}q;d" $dir/05_AnalysisOutput/beagle/pop_scaff.tsv`
pop=$(echo $ar | awk '{print $1}')
scaff=$(echo $ar | awk '{print $2}')

module load bcftools/1.9
bcftools view -r $scaff -S $dir/Sample_Lists/${pop}.ind $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.vcf.gz | \
  perl -pe "s/\s\.:/\t.\/.:/g" | \
  bcftools view -O z -o $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.${pop}_${scaff}.vcf.gz

module load beagle/5.4
beagle gt=$dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.${pop}_${scaff}.vcf.gz \
       out=$dir/05_AnalysisOutput/beagle/${pop}_${scaff} \
       impute=false \
       window=0.5 \
       overlap=0.025 \
       iterations=40 \
       phase-states=1000 \
       nthreads=16
rm $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.${pop}_${scaff}.vcf.gz

module load vcftools
bcftools view $ind_$scaff.polarized.vcf.gz | \
  bcftools +fill-tags | \
  vcf-annotate -a $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.bed.gz \
               -d key=INFO,ID=AA,Number=1,Type=String,Description='Ancestral Allele' \
               -c CHROM,-,FROM,INFO/AA | \
    bcftools view -O z -o $ind_$scaff.polarized.AA.phased.vcf.gz
rm $pop_$scaff.polarized.vcf.gz
```
## Caclulate iHS, Rsb, & XP-EHH
### Calculate iHH values for individual scaffolds for each population - **REHH_iHH.R**
```R
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

.libPaths("/share/cdfwwildlife/Capel_Dedicated/R_lib/4.3.1")
library(rehh)
library(stringr)

# Read in data and generate a master list of dataframes containing per-pop iHH, iES, and inES for given scaff
setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/beagle")
file <- list.files(path = ".", pattern = paste0("^",args[2],"_",args[1],".polarized.AA.phased.vcf.gz$"))
chr_num <- str_extract(string = args[1], pattern = "[0-9]+")
hh <- data2haplohh(hap_file = file,
                   polarize_vcf = T,
                   vcf_reader = "data.table",
                   min_perc_geno.mrk = 80)
scan.chrm <- scan_hh(hh,
                     maxgap = 1056,  # calculated as mean+5SD of distribution of all SNP distances
                     discard_integration_at_border = T)
scan.chrm$CHR <- chr_num

# Write output to file
setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh/iHH_per_pop_scaff")
write.table(scan.chrm, file = paste0(args[2],"_",args[1],".scan.chrm.tsv"),
            quote = F, sep = "\t", col.names = T, row.names = F)
```

### Combine iHH across scaffolds for each population and calculate iHS - **REHH_iHS.R**
```R
#!/usr/bin/env Rscript
.libPaths("/share/cdfwwildlife/Capel_Dedicated/R_lib/4.3.1")

library(rehh)
library(stringr)

pops <- c("NYPRE", "NYPOST", "PAPRE", "PAPOST", "PRE", "POST", "NY", "PA")
scan.chrms <- list()
wgscan.iHH <- list()
wgscan.iHH.pops <- list()

setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh/iHH_per_pop_scaff")

# Read in iHH data
cat("Reading in data...\n")
for (i in 1:length(pops)) {
  print(pops[i])
  files <- list.files(".", pattern = paste0("^",pops[i]))
  for (j in 1:length(files)) {
  scaf <- sub(".scan.chrm.tsv", "", sub(paste0(pops[i],"_"), "", files[j]))
  print(scaf)
    scan.chrms[[j]] <- read.csv(files[j], header = T, sep = "\t")
    if (j == 1) {
      wgscan.iHH <- scan.chrms[[j]]
    } else {
      wgscan.iHH <- rbind(wgscan.iHH, scan.chrms[[j]])
    }
  }
  wgscan.iHH.pops[[i]] <- wgscan.iHH
}

# Calculate individual population iHS
setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh/iHS_per_pop")
wgscan.iHS.pops <- list()
cr.iHS.pops <- list()
cat("Calculating iHS...\n")
for (i in 1:length(wgscan.iHH.pops)) {
  print(pops[i])
  wgscan.iHS.pops[[i]] <- ihh2ihs(wgscan.iHH.pops[[i]],
                                  freqbin = 0.01, # or 92
                                  min_maf = 0.03333333)
  write.table(wgscan.iHS.pops[[i]], file = paste(pops[i],"_wgscan_iHS.tsv",sep=""), quote = F,
              sep = "\t", col.names = T, row.names = F)
  cr.iHS.pops[[i]] <- calc_candidate_regions(wgscan.iHS.pops[[i]],
                                             threshold = 3,
                                             pval = T,
                                             window_size = 10000,
                                             overlap = 1000,
                                             min_n_extr_mrk = 2)
  write.table(cr.iHS.pops[[i]], file = paste(pops[i],"_cr_iHS.tsv",sep=""), quote = F,
              sep = "\t", col.names = T, row.names = F)
}

cat("Generating iHS plots...\n")
chr.names <- unique(as.numeric(cr.iHS.pops[[1]]$CHR))
pdf("iHS_diagnostic_plots.pdf")
for (i in 1:length(wgscan.iHS.pops)) {
  print(freqbinplot(wgscan.iHS.pops[[i]],
              main = paste("uniHS w/i freq. bins", pops[i])))
  print(distribplot(wgscan.iHS.pops[[i]]$ihs$IHS,
              xlab = "iHS",
              main = paste("Genome-wide distribution", pops[i])))
  print(distribplot(wgscan.iHS.pops[[i]]$ihs$IHS,
              qqplot = T,
              xlab = "iHS"))
  print(manhattanplot(wgscan.iHS.pops[[i]],
                pval = T,
                threshold = 3,  # either 3 or 4* to be comparable to 5SD
                chr.name = chr.names,
                cr = cr.iHS.pops[[i]],
                main = paste("iHS w/p-value cutoff & candidates", pops[i])))
}
dev.off()
```

### Combine iHH across scaffolds for each population and calculate Rsb & XP-EHH - **REHH_Rsb.R**
```R
#!/usr/bin/env Rscript
.libPaths("/share/cdfwwildlife/Capel_Dedicated/R_lib/4.3.1")

library(rehh)
library(stringr)

setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh")
pops <- c("NYPRE", "NYPOST", "PAPRE", "PAPOST", "PRE", "POST", "NY", "PA")
pop.comps <- read.csv("pop_comps.tsv", header = F, sep = "\t")
pop.comps.n <- data.frame(pop1 = c(1,1,3,3,5,7), pop2 = c(2,4,2,4,6,8))
scan.chrms <- list()
wgscan.iHH <- list()
scan.pops <- list()

# Read in iHH data
setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh/iHH_per_pop_scaff")
cat("Reading in data...\n")
for (i in 1:length(pops)) {
  print(pops[i])
  files <- list.files(".", pattern = paste0("^",pops[i]))
  for (j in 1:length(files)) {
    scaf <- sub(".scan.chrm.tsv", "", sub(paste0(pops[i],"_"), "", files[j]))
    print(scaf)
    scan.chrms[[j]] <- read.csv(files[j], header = T, sep = "\t")
  }
  scan.pops[[i]] <- scan.chrms
}

# Calculate pairwise Rsb
setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh/Rsb_per_comp")
cat("Calculating pairwise Rsb...\n")
for (i in 1:nrow(pop.comps)) {
  print(paste(pop.comps[i,1], pop.comps[i,2]))
  for (j in 1:length(scan.pops[[i]])) {  # iterate throug scaffs
    print(paste("SUPER__",scan.pops[[i]][[j]][1,1],sep=""))
    Rsb.chrom.left <- ines2rsb(scan_pop1 = scan.pops[[pop.comps.n[i,1]]][[j]],
                               scan_pop2 = scan.pops[[pop.comps.n[i,2]]][[j]],
                               popname1 = pop.comps[i,1],
                               popname2 = pop.comps[i,2],
                               p.side = "left")
    if (j == 1) {
      wgscan.Rsb.left <- Rsb.chrom.left
    } else {
      wgscan.Rsb.left <- rbind(wgscan.Rsb.left, Rsb.chrom.left)
    }
  }
  cat("  Writing file...\n")
  write.table(wgscan.Rsb.left, file = paste("Rsb_",pop.comps[i,1],"-",pop.comps[i,2],".tsv", sep = ""),
   	       quote = F, sep = "\t", col.names = T, row.names = F)

  cat("  Calculating candidate regions...\n")
  cr.Rsb.left <- calc_candidate_regions(wgscan.Rsb.left,
				threshold = -100,
				window_size = 10000,
				overlap = 1000,
				join_neighbors = F,
				min_n_extr_mrk = 0)
  cr.Rsb.left$Pop1 <- pop.comps[i,1]
  cr.Rsb.left$Pop2 <- pop.comps[i,2]
  cat("  Writing file...\n")
  write.table(cr.Rsb.left, file = paste("Rsb_10Kbwin_",pop.comps[i,1],"-",pop.comps[i,2],".tsv", sep = ""),
   	      quote = F, sep = "\t", col.names = T, row.names = F)
}

# Calculate pairwise XP-EHH
setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh/XP-EHH_per_comp")
cat("Calculating pairwise XP-EHH...\n")
for (i in 1:nrow(pop.comps)) {
  print(paste(pop.comps[i,1], pop.comps[i,2]))
  for (j in 1:length(scan.pops[[i]])) {  # iterate throug scaffs
    print(paste("SUPER__",scan.pops[[i]][[j]][1,1],sep=""))
    xpehh.chrom.left <- ies2xpehh(scan_pop1 = scan.pops[[pop.comps.n[i,1]]][[j]],
                               scan_pop2 = scan.pops[[pop.comps.n[i,2]]][[j]],
                               popname1 = pop.comps[i,1],
                               popname2 = pop.comps[i,2],
                               p.side = "left")
    if (j == 1) {
      wgscan.xpehh.left <- xpehh.chrom.left
    } else {
      wgscan.xpehh.left <- rbind(wgscan.xpehh.left, xpehh.chrom.left)
    }
  }
  cat("  Writing file...\n")
  write.table(wgscan.xpehh.left, file = paste("XP-EHH_",pop.comps[i,1],"-",pop.comps[i,2],".tsv", sep = ""),
   	      quote = F, sep = "\t", col.names = T, row.names = F)

  cat("  Calculating candidate regions...\n")
  cr.xpehh.left <- calc_candidate_regions(wgscan.xpehh.left,
					  threshold = -100,
					  window_size = 10000,
					  overlap = 1000,
					  join_neighbors = F,
					  min_n_extr_mrk = 0)
  cr.xpehh.left$Pop1 <- pop.comps[i,1]
  cr.xpehh.left$Pop2 <- pop.comps[i,2]
  cat("  Writing file...\n")
  write.table(cr.xpehh.left, file = paste("XP-EHH_10Kbwin_",pop.comps[i,1],"-",pop.comps[i,2],".tsv", sep = ""),
   	      quote = F, sep = "\t", col.names = T, row.names = F)
}
```

### Subset smoothed Rsb & XP-CLR down to 10% overlap (90% overlap generated) - **REHH_Rsb_filt.sbatch**
```shell
cd /share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh/Rsb_per_comp
for file in `ls *10Kb*`
do
    echo -e "Filtering $file ..."
    pre=$(echo $file | sed -E 's/.tsv//')
    head -n 1 $file > ${pre}_10perc_overlap.tsv
    tail -n +2 $file | cut -f 1 | uniq | sort -g | while read chrom
    do
        cat $file | awk -v chrom=$chrom '{if ($1 == chrom) print $0}' | head -n 1
        cat $file | awk -v chrom=$chrom '{if ($1 == chrom) print $0}' | tail -n +2 | awk 'NR % 9 == 0'
    done >> ${pre}_10perc_overlap.tsv
done

cd ../XP-EHH_per_comp
for file in `ls *10Kb*`
do
    echo -e "Filtering $file ..."
    pre=$(echo $file | sed -E 's/.tsv//')
    head -n 1 $file > ${pre}_10perc_overlap.tsv
    tail -n +2 $file | cut -f 1 | uniq | sort -g | while read chrom
    do
        cat $file | awk -v chrom=$chrom '{if ($1 == chrom) print $0}' | head -n 1
        cat $file | awk -v chrom=$chrom '{if ($1 == chrom) print $0}' | tail -n +2 | awk 'NR % 9 == 0'
    done >> ${pre}_10perc_overlap.tsv
done
```
