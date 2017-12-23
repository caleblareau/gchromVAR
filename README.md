# gchromVAR

### About:

Two outstanding challenges in the post-GWAS era are (1) the precise identification of causal variants within associated loci and (2) determining the exact mechanisms by which these variants result in the observed phenotypes, starting with identification of the pertinent cell types. To address (1), we used robust genetic fine mapping to identify hundreds of likely causal variants for 16 blood cell traits, allowing for up to 5 causal variants in each locus. We combined our fine-mapped results with high resolution open chromatin data for 18 primary hematopoietic populations and derived functional annotations to identify predicted target genes, mechanisms, and disease relevance. Moreover, we elucidate compelling anecdotes for the utility of this approach. To address (2), we developed a novel enrichment method (**gchromVAR**) that can discriminate between closely related cell types and score single cells for GWAS enrichment. 

We've implemented **gchromVAR** as an `R` package for computing cell-type specific GWAS enrichments from GWAS summary statistics and quantitative epigenomic data. This web resource and vignette compiliation shows how to reproduce these results in hematopoesis and how to run **gchromVAR** on other data sets. 

### Installation:

Once all of the dependencies for `gchromVAR` are installed, the package can be installed 
directly from GitHub by typing the following into an `R` console:

```
devtools::install_github("caleblareau/gchromVAR")
```


### Application to hematopoesis

We performed our genome-wide association studies (GWASs) on the UK Biobank (UKBB), which consists of 113,000 individuals of predominantly European descent in which 16 blood cell traits have been measured. These blood cell traits represent several important and distinct hematopoietic lineages, including red blood cells (RBCs), platelets, lymphocytes (T and B cells), monocytes, eosinophils, and basophils (bottom Figure 1). We overlaid these GWAS results with chromatin accessibility data derived from ATAC-seq (shown in the population tree). 

<p align="center">
  <br><br><br>
  <img src="vignettes/media/overview_Heme.png" width="40%"/><br>
  <b>Figure 1</b>. Overview of hematopoesis with cell types (top) and GWAS traits (bottom) explored thus far with <b>gchromVAR</b>. 
</p><br><br>

Although several methods have been developed to calculate enrichment of genetic variation with genomic annotations such as changes in chromatin accessibility (Trynka et al. 2015, Finucane et al. 2015), a method which takes into account both (1) the strength and specificity of the genomic annotation and (2) the probability of variant causality, accounting for LD structure, is needed to resolve associations within the stepwise hierarchies that define hematopoiesis. To these ends, we developed a new approach called genetic-chromVAR (**gchromVAR**), an adaptation of a recently described method, to measure the enrichment of regulatory variants in each cell state using our fine-mapped genetic variants and quantitative genomic annotations (Fig. 2A). Briefly, this method weights chromatin features by variant posterior probabilities and computes the enrichment for each cell type versus an empirical background matched for GC content and feature intensity. We show that g-chromVAR successfully identifies true enrichments of causal variants and is generally robust to variant posterior probability thresholds.

<br>
<p align="center">
  <br>
  <img src="vignettes/media/overview_Fig1.png" width="80%"/><br>
  <img src="vignettes/media/overview_Fig2.png" width="80%"/><br>
  <b>Figure 2</b>. Schematic and overview of results when applying <b>gchromVAR</b> to bulk populations. 
</p>
    
We applied **gchromVAR** to each of the 16 traits and 18 bulk ATAC-seq hematopoietic progenitor populations primarily sorted from the bone marrow of multiple healthy individual donors (Fig. 1A). We compared g-chromVAR to two state of the art methods: LDSR (Finucane et al. 2015), which calculates the enrichment for genome-wide heritability using binary annotations after accounting for LD and overlapping annotations, and goShifter (Trynka et al. 2015), which calculates the enrichment of tight LD blocks containing sentinel GWAS SNVs for binary annotations. Using a Bonferroni correction, g-chromVAR identified 22 trait-tissue enrichments, LDSR identified 71, goShifter identified 39, and chromVAR identified 79 (Fig. 2C).

In order to compare the performance of enrichment tools, we leveraged our knowledge of the hematopoietic system and devised a lineage specificity test. For any measured cell trait we identified all possible upstream progenitors that could be passed through before terminal differentiation (Fig. 1A). For example, the differentiation of an RBC is generally thought to begin at the hematopoietic stem cell (HSC) and progress through multipotent progenitor (MPP), common myeloid progenitor (CMP), and megakaryocyte erythroid progenitor (MEP) before reaching the erythroid progenitor (Ery) stage. Thus, the lineage specificity test is a nonparametric rank-sum test that measures the relative ranking of lineage specific trait-cell type pairs relative to the non lineage specific traits for each of the compared methodologies. Using this metric for specificity, we found that g-chromVAR vastly outperformed all three other methods (Fig 2D). Additionally, we found that 21/22 (95%) of g-chromVAR trait-cell type enrichments were supported by LDSR, all of which were lineage specific (Fig. 2C). For certain traits such as monocyte count, we found highly similar enrichment patterns for g-chromVAR and LDSR, but non lineage enrichments for **chromVAR**. For other traits, such as mean reticulocyte volume, g-chromVAR identified only the two most terminally proximal cell types (MEP and Ery) as significantly enriched for the trait, whereas LDSR non-specifically identified 15/18 of the investigated cell types as enriched after Bonferroni correction. We note that we can improve the lineage specificity of LDSR by including all hematopoietic ATAC-seq annotations in the model as covariates, but this results in a loss of power.

Having validated our approach, we investigated cell type enrichments for each of the 16 traits. We found that the most lineage-restricted progenitor populations were typically most strongly enriched (Fig. 2E-H). For example, RBC count was most strongly enriched in erythroid progenitors (Fig. 3E), platelet count was most strongly enriched in megakaryocytes (platelet progenitors) (Fig. 2F), and lymphocyte count was most strongly enriched in CD4+ and CD8+ T cells (Fig. 2H). In several instances, we observed significant enrichments for traits in earlier progenitor cell types within each lineage, including enrichment for platelet traits in CMPs and enrichment for monocyte traits in a specific subpopulation of GMPs. Building on several studies that recently demonstrated transcriptomic and chromatin accessibility heterogeneity within these populations, we next sought to apply g-chromVAR to single-cell ATAC-seq (scATAC-seq) data in order to interrogate the impact of common genetic variation underlying blood cell traits in heterogeneous cell populations.


### Contact:
[Caleb Lareau](mailto:caleblareau@g.harvard.edu) developed and maintains this package with Erik Bao and Jacob Ulirsch. 

<br><br>