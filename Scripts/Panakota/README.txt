##### Bacteroides de novo protein prediction, phylogenetic tree, orthologs search, LexA copy number identification

# downloaded 51 representative Bacteroides genomes from refseq (fna)
# PANACOTA
# installation
singularity pull --name panacota.img docker://gempasteur/panacota
# PanACoTA - v. 1.3.1

/scratch/okhtienk/panacota.img prepare --norefseq -o Bacteroides_51_prep -d /home/okhtienk/Bacteroides_51 --max_dist 1.0 --l90 105
# --min_dist MIN_DIST   By default, genomes whose distance to the reference is not between 1e-4 and 0.06 are discarded. You can specify your own lower limit (instead of 1e-4) with this option.
# --max_dist MAX_DIST   By default, genomes whose distance to the reference is not between 1e-4 and 0.06 are discarded. You can specify your own lower limit (instead of 0.06) with this option.
# L90 < 100, n_contigs < 999 - default
/scratch/okhtienk/panacota.img annotate --info ~/Bacteroides_51_prep/LSTINFO-NA-filtered-0.0001_1.0.txt -r ./Panakota_50 -n BACT --threads 20


# Got 50 assembly (1 genome was discarded)
# renamed files using python code from Panakota_output_rename_for_tree.ipynb

# OrthoFinder version 2.5.4
python OrthoFinder_source/orthofinder.py -f Panakota_proteins_50

# Miltiple alignment of orthologous proteins from previous launch of OrthoFinder (Bacteroides only, annotation from ncbi) with Clustal Omega online. Save in STOCKHOLM format.
# Launch HMMER-3.3.2

../src/hmmbuild bact_only.hmm ~/Bacteroides_uniformis_annotations/HMM/aln-stockholm.txt

# searching in B.theta protein multi-fasta (to find newly assigned id)

../src/hmmsearch bact_only.hmm ~/Bacteroides_uniformis_annotations/HMM/Bacteroides_thetaiotaomicron_GCF_014131755.1.faa > lexA_hits.out
