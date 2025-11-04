#!/bin/bash

echo "Copy and run commands directly on command line in interactive session!!"

module load minimap2
module load samtools

# Align genomes allowing 0.1% divergence
minimap2 -x asm5 $ref $out > ${out}/${pre}.paf  # create file for plotting
minimap2 -ax asm5 $ref $out | samtools view -bh -@ 2 -m 3G - | samtools sort -@ 2 -m 3G -O bam > ${out}/${pre}.bam
