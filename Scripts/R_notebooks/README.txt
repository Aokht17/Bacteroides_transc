R Rmd notebooks description:

- WGCNA: code for WGCNA analysis of transcriptomics and proteoimics data

DGE folder

- allmice_VS_allcult.Rmd: transcriptomics, B.uniformis data, DGE all samples from mice VS samples from culture (without treatment)
- antib_deseq.Rmd: transcriptomics, B.uniformis data, DGE metronidazole VS DMSO treatment in culture
- b_theta_deseq.Rmd: transcriptomics, B.thetaiotaomicron data, DGE of wild type and LexA knockouts (tried linear models)
- Custom_GO_GSEA.Rmd: code to run GO enrichment with custom GO terms annotation table and adjustment to many un-annotated genes. Arranged in function in Modules_enrichment.Rmd
- exp_vs_mono.Rmd: transcriptomics, B.uniformis data, DGE monocolonized mice VS exponential growth phase in culture
- lag_VS_exp.Rmd: transcriptomics, B.uniformis data, DGE early exponential VS exponential growth phase
- Modules_enrichment.Rmd: contains custom function for GO analysis for sparsely annotated genomes and the 2nd for GO analysis using list of genes without LFC values
- Mono_vs_community_cecum.Rmd: transcriptomics, B.uniformis data, mice cecum content, DGE monocolonized mice VS 8-species community
- proteomics_DIF.Rmd: DGE analysis of B.uniformis proteomics data (RAW+DESeq2, not VSN+limma) and scatterplots comparison with transcriptomics results
- proteom_transc_scatterplots.Rmd: Scatterplots proteomics VS transcriptomics (VST values)
- stat_VS_exp.Rmd: transcriptomics, B.uniformis data, DGE stationary VS exponential growth phase
- WT_VS_bothKO_theta.Rmd: transcriptomics, B.thetaiotaomicron data, DGE wild type VS double LexA knockout
