# Code for Eastern *Myotis lucifugus* WNS Selection Manuscript
## Manuscript Information
**Title:** A multifaceted approach reveals complex genomic mediation of white-nose syndrome adaptative response in the little brown bat (*Myotis lucifugus*)

**Dataset author:** Samantha L. R. Capel; email: <Samantha.Capel@wildlife.ca.gov> or <slr.capel2@gmail.com>

**Citation:** Samantha L. R. Capel, Devaughn L. Fraser, Kenneth A. Field, DeeAnn M. Reeder, Amy L. Russell, Peter H. Sudmant, Juan Manuel Vazquez, Maarten J. Vonhof, Thomas M. Lilley, and Michael R. Buchalski. (2026). A multifaceted approach reveals complex genomic mediation of white-nose syndrome adaptative response in the little brown bat (*Myotis lucifugus*). *Molecular Ecology* **VOL**(ISS):PG-PG. https://doi.org/

## Data availability
Raw data: [NCBI BioProject PRJNA1353610](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1353610/), BioSample accession nos. SAMN52933055–SAMN52933113

Output data available through Dryad: https://doi.org/10.5061/dryad.ncjsxkt66

## Directories & Files
- **Call_SNPs_Pipeline** - scripts for calling & filtering SNPs; see the [dedicated GitHub repository](https://github.com/slcapel/DRAGEN-GATK4_SNP_calling_pipeline) for more detailed information
- **data_files** - metadata & auxiliary files
- ANGSD_Fst.md - calculate smoothed Fst using ANGSD
- ANGSD_Fst_outliers.R - call Fst outliers & generate Manhattan plots
- XP-CLR.md - run XP-CLR
- XP-CLR_outliers.R - call XP-CLR outliers & generate Manhattan plots
- rehh.md - scripts for running REHH & calculating Rsb
- REHH_Rsb_outliers.R - call Rsb outliers & generate Manhattan plots
- SNPRelate.R - calculate pairwise kinship on all individuals
