# A multifaceted approach reveals complex genomic mediation of white-nose syndrome  resistance  in the little brown bat (*Myotis lucifugus*)
Dataset DOI: [10.5061/dryad.ncjsxkt66](https://doi.org/10.5061/dryad.ncjsxkt66)

## Description of the data and file structure
This dataset contains all analytical code for this manuscript along with associated data files, two versions of the final VCF, and selection statistic output files. The genome annotation used with this data is available thorugh https://github.com/docmanny/myotis-gene-annotations.

## Data Files
#### VCFs:
- gatk.snp.qual_hard_filtered_autosomes.vcf.gz -- full set of filtered SNPs
- gatk.snp.qual_hard_filtered_autosomes_thin.vcf.gz -- filtered SNPs thinned by 10 Kbp distance

#### SNP Calling Pipeline Input Files:
- RG_info.tsv -- read group information for variant calling
- all_samples.txt -- sample IDs for all individual FASTQ files
- all_samples.merge.txt -- sample IDs merging all files for each individual

#### Selection Statistic Output Files:
- ANGSD_Fst_10kbp_windows.zip -- final output files from **ANGSD_Fst.md**
- XP-CLR_10kbp_windows.zip -- final output files from **XP-CLR_outliers.R**
- REHH_Rsb_10_Kbp_windows.zip -- final output files from **REHH_Rsb_outliers.R**

#### Sample and Population Text Files for Scripts:
- pops.txt -- text file of all population designations
- Population_Map.tsv -- population and sampling site information for each individual
- pop_comps.tsv -- text file of all pairwise comparisons
- PRE.ind -- all pre-WNS individuals
- pop_scaff.tsv -- text file input for phasing genotypes

#### Alignment and Variant Statistics:
- hicov_coverage.tsv -- per-scaffold alignment statistics for "high coverage" individuals
- lowcov_coverage.tsv -- per-scaffold alignment statistics for "shallow coverage" individuals
- SNP_filtering.txt -- Type and number of variants filtered per scaffold

## Scripts & Code
#### SNP Calling Pipeline:
- **call_SNPs_pipeline.zip**
- Components:
  - clean_align_callSNPs.sbatch -- master script
  - HTS_preproc.slurm -- clean FASTQs; pipeline component
  - hashDRAGMAP.slurm -- build reference genome hash table; pipeline component
  - alignDRAGMAP.slurm -- align FASTQs to reference genome; pipeline component
  - samtools_merge.slurm -- merge sample bams; pipeline component
  - genome_wins.slurm -- create .bed for genome windows; pipeline component
  - align_stats.slurm -- calculate alignment statistics; pipeline component
  - STRtable.slurm -- create STR table; pipeline component
  - bam_to_gvcf.slurm -- call variants per individual; pipeline component
  - gvcf_to_vcf_scaff.slurm -- joint call variants; pipeline component
  - vcf_scaff_to_snp.vcf.slurm -- combine VCFs and filter; pipeline component
- Dependencies: [HTStream](https://s4hts.github.io/HTStream/#hts_QWindowTrim), [DRAGMAP v1.2.1](https://github.com/Illumina/DRAGMAP), [picard](https://github.com/broadinstitute/picard), [GATK4](https://github.com/broadinstitute/gatk), [samtools](https://github.com/samtools/samtools), [bedtools](https://github.com/arq5x/bedtools2), and [bcftools](https://github.com/samtools/bcftools)
- Inputs: [raw sequencing files](SRA LINK), [*M. lucifugus* reference genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_048340685.1/), RG_info.tsv, all_samples.txt, all_samples.merge.txt
- For more details see the associated [github repository](https://github.com/slcapel/DRAGEN-GATK4_SNP_calling_pipeline)

#### Selection Statistic Calculation:
- **ANGSD_Fst.md** -- calculate population pariwise F<sub>ST</sub> in 10 Kbp windows for all population comparisons
  - Dependencies: [ANGSD v0.940](https://anaconda.org/bioconda/angsd) and [bcftools](https://github.com/samtools/bcftools)
  - Inputs: gatk.snp.qual_hard_filtered_autosomes.vcf.gz, pops.txt, Population_Map.txt, pop_comps.tsv
- **XP-CLR.sbatch** -- calculate XP-CLR in 10 Kbp windows for all population comparisons
  - Dependencies: [XP-CLR](https://github.com/hardingnj/xpclr), [bcftools](https://github.com/samtools/bcftools), and [tabix](https://www.htslib.org/doc/tabix.html)
  - Inputs: gatk.snp.qual_hard_filtered_autosomes.vcf.gz, pop_comps.tsv
- **rehh.md** -- polarize VCF alleles, phase genotypes, and calculate Rsb 10 Kbp windows for all population comparisons
  - Dependencies: [java 17](https://anaconda.org/conda-forge/openjdk/files?page=0&type=&version=17.0.3), [bcftools](https://github.com/samtools/bcftools), [vcftools](https://vcftools.sourceforge.net/man_latest.html), [vcffilterjdk](https://lindenb.github.io/jvarkit/VcfFilterJdk.html), [beagle v5.4](https://faculty.washington.edu/browning/beagle/beagle_5.4_18Mar22.pdf), [R](https://www.r-project.org/), [rehh](https://cran.r-project.org/web/packages/rehh/index.html), and [strigr](https://cran.r-project.org/web/packages/stringr/index.html)
  - Inputs: gatk.snp.qual_hard_filtered_autosomes.vcf.gz, PRE.ind, pop_scaff.tsv, pop_comps.tsv

#### Identify Selection Statistic Outlier Windows:
- **ANGSD_Fst_outliers.R** -- identify outlier windos from final **ANGSD_Fst.md** output
  - Dependencies: [R](https://www.r-project.org/) and [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
  - Input: ANGSD_Fst_10kbp_windows.zip (unzipped) 
- **XP-CLR_outliers.R** -- identify outlier windos from final **XP-CLR.sbatch** output
  - Dependencies: [R](https://www.r-project.org/) and [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
  - Input: XP-CLR.zip (unzipped)
- **REHH_Rsb_outliers.R** -- identify outlier windos from final **rehh.sbatch** output
  - Dependencies: [R](https://www.r-project.org/) and [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
  - Input: REHH_Rsb_10_Kbp_windows.zip (unzipped)

### Other Scripts:
- **SNPRelate.R** -- Calculate pairwise kinship using SNPRelate's IBD MLE method
  - Requirements: [R](https://www.r-project.org/), [SNPRelate](https://www.rdocumentation.org/packages/SNPRelate/versions/1.6.4), and [gdsfmt](https://www.bioconductor.org/packages/release/bioc/html/gdsfmt.html)
  - Input: gatk.snp.qual_hard_filtered_autosomes_thin.vcf.gz
