######################################

Dir content description:

Panakota
R_notebooks
jupyter_notebooks

- annogesic.txt: pipeline for dRNA-seq analysis
- egg_gene_names.tsv & index_matching.csv: tables with genes names and IDs (my annotation, modified EBI+EggNOG+UNIPROT)
- meta_salmon.sh: transcriptomics analysis using alignment on transcriptome and Salmon tool (fast)
- metatranscripromics_pipe.sh: long pipeline using genome and hisat2 (slow)
- Motifs_regNet.py: summary of naive functions that can be used for binding motifs investigation. See jupyter notebooks for more
- N0.tsv: root level of OrthoFinder prediction results. 'LexA' orthologs found in different Bacteroides species
- phylogenetic_tree_Bacteroides.pdf: picture of tree. Red dots show the number of found LexA gene copies. The tree was built on the single copy orthologous genes by the OrthoFinder tool. The tree is rooted by STRIDE. Visualized in iTOL 
- prepDE.py: script (https://github.com/gpertea/stringtie) needed for the metatranscripromics_pipe.sh launch


