```bash
#!/bin/bash

aklog
source /share/cdfwwildlife/Capel_Dedicated/miniconda/etc/profile.d/conda.sh
conda activate XP-CLR
module load bcftools

# NOTE: index vcf using bcftools prior to executing 

dir="/share/cdfwwildlife/MYLU_NovaSeq"  # root directory
vcf="${dir}/04_GATKvcfs/gatk.snp.super_filtered_autosomes_thin.vcf.gz"  # vcf to be analyzed
pop1="${dir}/Sample_Lists/NYPRE.ind"  # text file with one sample name per line
pop2="${dir}/Sample_Lists/NYPOST.ind"
pop3="${dir}/Sample_Lists/PAPRE.ind"
pop4="${dir}/Sample_Lists/PAPOST.ind"
outdir="${dir}/05_AnalysisOutput/XP-CLR"

c=0
ar1=($pop1 $pop2 $pop3)  # the first n-1 population list variables
ar2=($pop2 $pop3 $pop4)  # the last n-1 population list variables

tabix $vcf

for i in {0..2}
do
    for j in $(seq $c 2)
    do
        bcftools index -s $vcf | cut -f 1 | while read scaff
        do
            p1=${ar1[$i]}
            p1=$(echo ${p1##*/})
            p1=$(echo ${p1%.*})
            p2=${ar2[$i]}
            p2=$(echo ${p2##*/})
            p2=$(echo ${p2%.*})
            echo "RUNNING XP-CLR ON $scaff for $p1 vs $p2..."
            xpclr -F vcf -I $vcf -C $scaff -Sa ${ar1[$i]} -Sb ${ar2[$j]} --size 50000 --step 45000 --minsnps 2 -O $outdir/${scaff}_${p1}_${p2}.out -V 10
        done
    done
    c=$((c+1))
done#!/bin/bash

#SBATCH --job-name=xpclr
#SBATCH --time=30-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --array=1-168
#SBATCH --output=/share/cdfwwildlife/MYLU_NovaSeq/SCRIPTS/slurmout/xpclr/%x_%A_%a.out
#SBATCH --error=/share/cdfwwildlife/MYLU_NovaSeq/SCRIPTS/slurmout/xpclr/%x_%A_%a.err

dir="/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/XP-CLR"
comp=`sed "${SLURM_ARRAY_TASK_ID}q;d" $dir/comps.tsv`
scaff=$(echo ${comp} | awk '{print $1}')
pop1=$(echo ${comp} | awk '{print $2}')
pop2=$(echo ${comp} | awk '{print $3}')
p1=$(echo ${comp} | awk '{print $4}')
p2=$(echo ${comp} | awk '{print $5}')
start=`date +%s`
echo $HOSTNAME
echo -e "${scaff}\t${p1}\t${p2}"

#initiate conda for XP-CLR
aklog
source /share/cdfwwildlife/Capel_Dedicated/miniconda/etc/profile.d/conda.sh
conda activate XP-CLR
aklog
module load bcftools

# NOTE: index vcf using bcftools & tabix prior to executing

dir="/share/cdfwwildlife/MYLU_NovaSeq"  # root directory
vcf="${dir}/gatk.snp.qual_hard_filtered_autosomes.vcf.gz"  # vcf to be analyzed
outdir="${dir}/10Kb_500SNPsmax"
[[ -d $outdir ]] || mkdir -p $outdir

xpclr -F vcf -I $vcf -C $scaff -Sa $pop1 -Sb $pop2 --size 10000 --step 9000 --minsnps 2 --maxsnps 500 -O $outdir/${scaff}_${p1}_${p2}.out -V 10

(exit) && echo success
```
