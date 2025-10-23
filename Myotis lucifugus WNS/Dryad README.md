# A multifaceted approach reveals complex genomic mediation of white-nose syndrome  resistance  in the little brown bat (*Myotis lucifugus*)
Dataset DOI: [10.5061/dryad.ncjsxkt66](https://doi.org/10.5061/dryad.ncjsxkt66)

## Description of the data and file structure
This dataset contains all analytical code for this project along with associated files, two versions of the final VCF, and selection statistic output files.

## Data Files
#### VCFs:
- gatk.snp.qual_hard_filtered_autosomes_thin.vcf.gz
- gatk.snp.qual_hard_filtered_autosomes.vcf.gz

#### SNP Calling Pipeline Input Files:
- RG_info.tsv
- all_samples.txt
- all_samples.merge.txt

#### Selection Statistic Output Files:
- ANGSD_Fst_10kbp_windows.zip
- XP-CLR.zip
- REHH_Rsb_10_Kbp_windows.zip

#### Other Files:
- SNP_filtering.txt -- Type and number of variants filtered per scaffold

## Scripts & Code
#### SNP Calling Pipeline:
- **call_SNPs_pipeline.zip**
- Components:
  - X
- Requirements: 
- Input(s): 

#### Selection Statistic Calculation:
- **ANGSD_Fst.md**
  - Requirements: ANGSD v0.940, bcftools, 
  - Input(s): gatk.snp.qual_hard_filtered_autosomes.vcf.gz, pops.txt, Population_Map.txt, pop_comps.tsv
- **XP-CLR.sbatch**
  - Requirements: 
  - Input(s): 
- **rehh.md**
  - Requirements: 
  - Input(s): 

#### Identify Selection Statistic Outliers:
- **ANGSD_Fst_outliers.R**
  - Requirements: 
  - Input(s): 
- **XP-CLR_outliers.R**
  - Requirements: 
  - Input(s): 
- **REHH_Rsb_outliers.R**
  - Requirements: 
  - Input(s): 

### Other Scripts:
- **SNPRelate.R** -- Calculate pairwise kinship using SNPRelate's IBD MLE method
  - Requirements: [R](https://www.r-project.org/), [SNPRelate](https://www.rdocumentation.org/packages/SNPRelate/versions/1.6.4), and [gdsfmt](https://www.bioconductor.org/packages/release/bioc/html/gdsfmt.html)
  - Input(s): gatk.snp.qual_hard_filtered_autosomes_thin.vcf.gz
