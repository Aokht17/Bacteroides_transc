#!/bin/bash

# needs >= 32 cores
# hint: you can create transcriptome file with gffread -F -w transcriptome.fa -g genome.fna annot.gff
raw=$1
ref_transc=$2

raw_files="$(find $raw -type f )"

############## MAIN ####################

mkdir raw_QC
cd raw_QC
fastqc $raw_files -o .
multiqc .
echo "QC is ready"
cd ..

mkdir Trim
cd Trim
mkdir trim_QC
trim_galore --fastqc_args "--outdir ./trim_QC" --length 50 --cores 4 --trim-n $raw_files
cd trim_QC
multiqc .
echo "Trimming is ready"
cd ../..

mkdir Salmon
cd Salmon
salmon index -t $ref_transc -i salmon_idx
trim_path=../Trim
for file in $(find $trim_path -mindepth 1 -type f -name '*trimmed.fq'); do
	samp=`basename ${file}`
	echo "Processing sample ${samp}"
	salmon quant -i salmon_idx -l A -r $file -p 32 --validateMappings -o quants/${samp}_quant
done 
cd ../


echo "COMPLETED"


