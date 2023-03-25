# Computation-Analysis-of-TCGA-Cancer-Data

The project involves computational analysis of large data sets containing biomolecular measurements annotated with phenotypes from numerous anonymized cancer patients. These data sets are among those that have now become publicly available from TCGA (The Cancer Genome Atlas) https://cancergenome.nih.gov. Such analysis provides unique opportunities for discovering associations leading to biological hypothesis, which, when validated, have the potential to lead to diagnostic, prognostic and therapeutic products.

Specifically, this project focuses on analyzing gene expression data to derive associations among genes as well as disease phenotypes, co-expression signatures, and survival analysis. I also perform literature search in the biological/medical literature to understand the underlying biology behind results reached by computational methods.

The MATLAB workspace containing publicly available TCGA data sets with various kinds of data from cancer biopsies (e.g., gene expression, methylation, somatic mutations, survival, and other phenotypes):
1. BRCA (Breast)
2. COAD (Colon)
3. GBM (Glioblastoma)
4. HNSC (Head & Neck)
5. KIRC (Kidney)
6. LAML (Leukemia)
7. LUAD (Lung Adenocarcinoma)
8. LUSC (Lung Squamous Cell) -> I focused on this cancer type in this project
9. OV (Ovarian)
10. UCEC (Endometrial)
In addition to gene expression before log normalization (variables end in “_ge”, it contains the DNA methylation values for 20,225 genomic locations from all patients in ten cancer types (methylation variables end in “_me”. The variable Probe contains the official names of the methylation probes (starting from “cg,” and the variable Gene_meth contains the corresponding names of genes at those locations. Note that a gene may be associated with multiple methylation probes, and there are methylation probes that are not associated with genes. As discussed when referring to CpG islands in the genome, DNA methylation is an epigenetic mechanism that can control gene expression in various ways, usually methylation tends to silence the corresponding gene by inhibiting its expression “switching” it to an “off” position.
