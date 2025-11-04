#!/bin/bash

#SBATCH --job-name=busco_comp
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --ntasks-per-node=48
#SBATCH --array=1-9
#SBATCH --output=/share/cdfwwildlife/Golden_Eagle_RG_CCGP/SCRIPTS/slurmout/%x_%A_%a.out
#SBATCH --error=/share/cdfwwildlife/Golden_Eagle_RG_CCGP/SCRIPTS/slurmout/%x_%A_%a.err

start=`date +%s`
echo $HOSTNAME

dir="/share/cdfwwildlife/Golden_Eagle_RG_CCGP"
ortho="aves_odb10"
genomes=($dir/genomes/*.fna)
genome=${genomes[$SLURM_ARRAY_TASK_ID-1]}
pre=$(basename ${genome} .fna)
[[ -d $dir/busco ]] || mkdir -p $dir/busco  # create directory for busco results for all genomes
[[ -d $dir/busco/${pre}_${ortho} ]] || mkdir -p $dir/busco/${pre}_${ortho}  # create subdirectories for each genome

cd $dir/genomes
aklog
source /share/cdfwwildlife/Capel_Dedicated/miniconda/etc/profile.d/conda.sh
conda activate busco5

echo ${genome}
echo "Running BUSCO analysis..."
busco --offline \
      -f \
      -c 48 \
      -m geno \
      -l $dir/${ortho} \
      -i ${genome} \
      -o ${pre}_${ortho} \
      --out_path $dir/busco

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime"
