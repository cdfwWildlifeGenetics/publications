## SNP Calling Pipeline Input Files:
- RG_info.tsv -- read group information for variant calling
- all_samples.txt -- sample IDs for all individual FASTQ files
- all_samples.merge.txt -- sample IDs merging all files for each individual

## Sample and Population Text Files for Scripts:
- pops.txt -- text file of all population designations
- Population_Map.tsv -- population and sampling site information for each individual
- pop_comps.tsv -- text file of all pairwise comparisons
- PRE.ind -- all pre-WNS individuals
- pop_scaff.tsv -- text file input for phasing genotypes

## Alignment and Variant Statistics:
- hicov_coverage.tsv -- per-scaffold alignment statistics for "high coverage" individuals
- lowcov_coverage.tsv -- per-scaffold alignment statistics for "shallow coverage" individuals
- SNP_filtering.txt -- Type and number of variants filtered per scaffold

## Large data files archived through [Dryad](https://doi.org/10.5061/dryad.ncjsxkt66)
- Full and thinned SNP dataset VCFs
- Fst, XP-CLR, and Rsb analysis output data
