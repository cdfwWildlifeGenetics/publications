#!/bin/bash

echo "Copy and run commands directly on command line in interactive session!!"

module load minimap2
module load samtools

# Align genomes allowing 0.1% divergence
minimap2 -x asm5 $ref $out > ${out}/AqChCa_to_AqChCh_align.paf  # create file for plotting
minimap2 -ax asm5 $ref $out | samtools view -bh -@ 2 -m 3G - | samtools sort -@ 2 -m 3G -O bam > ${out}/AqChCa_to_AqChCh_align.bam

# Filter alignment
cat ${out}/AqChCa_to_AqChCh_align.paf | \
  awk '{if ($12 > 40 && $11 > 50000 && $2 > 15000000 && $7 > 15000000) print $0}' | \
  awk '{if ($5 == "-") print $6" "$8" "$9" "$1" "$4" "$3; else print $6" "$8" "$9" "$1" "$3" "$4}' > AqChCa_to_AqChCh_align.txt
  
