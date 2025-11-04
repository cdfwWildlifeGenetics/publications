#!/bin/bash

#SBATCH --job-name=liftoff
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 10-00:00:00
#SBATCH --mem=80GB
#SBATCH --output=/share/cdfwwildlife/Golden_Eagle_RG_CCGP/SCRIPTS/slurmout/%x_%A.out
#SBATCH --error=/share/cdfwwildlife/Golden_Eagle_RG_CCGP/SCRIPTS/slurmout/%x_%A.err

start=`date +%s`
echo $HOSTNAME

con="/share/cdfwwildlife/Capel_Dedicated/miniconda"
dir="/share/cdfwwildlife/Golden_Eagle_RG_CCGP"

aklog
source ${con}/etc/profile.d/conda.sh
conda activate base

liftoff -g $dir/genomes/GCF_900496995.4_bAquChr1.4_genomic.gtf -polish  -chroms $dir/minimap2/chroms.txt -dir $dir/liftoff -o $dir/liftoff/AqChCh_to_AqChCa_liftover -p 16 $dir/genomes/GCA_030820615.1_bAquChr2.0.hap1_genomic.fna $dir/genomes/GCF_900496995.4_bAquChr1.4_genomic.fna

for i in `ls Aq*`
do
    mv $i $i.gtf
done

end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime"
