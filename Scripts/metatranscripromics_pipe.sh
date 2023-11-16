#!/bin/bash

# takes long, needs >= 32 cores
raw=$1
ref_fasta=$2
ref_gtf=$3

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

mkdir Align
cd Align
mkdir sorted
trim_path=../Trim
hisat2-build -p 32 $ref_fasta idx
for file in $(find $trim_path -mindepth 1 -type f -name '*trimmed.fq'); do
	echo $(basename $file)
	hisat2 -x idx -U $file -S $(basename $file).sam
	samtools view -S -b $(basename $file).sam | samtools sort > sorted/$(basename $file)_sorted.bam 
done
echo "alignment is done"
cd sorted
for file in $(find ./ -mindepth 1 -type f -name '*_sorted.bam'); do
	bamtools split -reference -in $file
	FILE=$(echo "$(basename $file)" | cut -d "_" -f 1)
	samtools merge -@ 16 ${FILE}_final.bam *REF_contig_*
	rm *_sorted.REF_*
done 
cd ../..


mkdir Count
cd Count
mkdir reest
align_path=../Align/sorted
for file in $(find $align_path -mindepth 1 -type f -name '*final.bam'); do
	stringtie -p 16 -G $ref_gtf -l $(basename $file) -o $(basename $file).gtf $file
done

ls *.gtf > mergelist.txt
stringtie --merge -p 8 -G $ref_gtf -o stringtie_merged.gtf mergelist.txt
gffcompare -r $ref_gtf -G -o merged stringtie_merged.gtf
echo "reestimation after merging"

for file in $(find $align_path -mindepth 1 -type f -name '*final.bam'); do
	stringtie -e -B -p 16 -G stringtie_merged.gtf -o reest/$(basename $file).gtf $file
done
cd reest
for file in $(find ./ -mindepth 1 -type f -name '*.gtf'); do
	FILE=$(echo "$(basename $file)" | cut -d "_" -f 1)
	printf "%s\t%s\n" "$FILE" "$(basename $file)" >> samples.txt
done
python3 /scratch/okhtienk/metatrans/prepDE.py -i samples.txt
cd ..
echo "COMPLETED"


