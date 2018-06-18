# Maternal group 2 innate lymphoid cells control fetal growth and protect from endotoxin-induced abortion in mice #

**Elisa Balmas<sup>1,2,3,‡</sup>, Batika M. J. Rana <sup>4,‡</sup>, Russell S. Hamilton<sup>2,3</sup>, Norman Shreeve<sup>1,3</sup>, Jens Kieckbusch<sup>1,3</sup>, Irving Aye<sup>1,3</sup>, Delia A. Hawkes<sup>1</sup>, Sophie Trotter<sup>1</sup>, Jorge López-Tello2,3, Hannah Yong<sup>2,3</sup>, Salvatore Valenti<sup>1</sup>, Amanda Sferruzi-Perri<sup>2,3</sup>, Francesca Gaccioli<sup>1,3</sup>, Andrew N. J. McKenzie<sup>4,‡</sup> and Francesco Colucci<sup>1,3,‡,§</sup>**

<sup>1</sup> Department of Obstetrics and Gynaecology, University of Cambridge, University of Cambridge School of Clinical Medicine, NIHR Cambridge Biomedical Research Centre <br>
<sup>2</sup> Department of Physiology, Development and Neuroscience, University of Cambridge<br>
<sup>3</sup> Centre for Trophoblast Research, University of Cambridge<br>
<sup>4</sup> MRC Laboratory of Molecular Biology, Francis Crick Avenue, Cambridge Biomedical Campus, Cambridge CB2 0QH, UK<br>
<sup>‡</sup> Equal contribution<br>
<sup>§</sup> Corresponding author: fc287@medschl.cam.ac.uk <br>

### Citation ###

Balmas, E., Rana, B.M.J., Hamilton, R.S., Shreeve, N., Kieckbusch, J., Aye, I., Hawkes, D.A., Trotter, S., López-Tello, J., Yong, H., Valenti., S., Sferruzi-Perri, A., Gaccioli, F., McKenzie, A.N.J. & Colucci, F. (2018) Maternal group 2 innate lymphoid cells control fetal growth and protect from endotoxin-induced abortion in mice. [bioRxiv 348755](https://doi.org/10.1101/348755)

### Abstract ###
Group 2 innate lymphoid cells (ILC2s) adapt to tissue physiology and contribute to immunity, inflammatory pathology and metabolism. We show that mouse uterine ILC2s have a heightened type-2 gene signature and expand during pregnancy. Indeed, maternal ILC2s protect against fetal growth restriction and mortality upon systemic endotoxin challenge. The absence of ILC2s leads to utero-placental abnormalities, including poor vascular remodelling, increased Il1b and decreased Il4, Il5, and Il13 gene expression, and reduced dendritic cells (DCs) and macrophage alternative activation. Placentas exhibit signs of adaptation to stress, including larger maternal blood spaces and increased expression of nutrient transporter genes. Endotoxin induces the expansion of IL-1-producing uterine DCs and, in response, more uterine ILC2s produce IL-4, IL-5 and IL-13.  In a protective feedback mechanism, these cytokines suppress IL-1-producing uterine DCs, in line with a protective role of uILC2s against endotoxin-induced abortion. Uterine ILC2s emerge as pivotal for both normal and complicated pregnancies.

### Data Processing ###
Data were aligned to GRCm38 mouse genome (Ensembl Release 86) with TopHat2 (v2.1.1, using bowtie2 v2.2.9) with a double map strategy. Alignments and QC were processed using custom ClusterFlow (v0.5dev) pipelines and assessed using MultiQC (0.9.dev0). Gene quantification was determined with HTSeq-Counts (v0.6.1p1). Additional quality control was performed with feature counts (v 1.5.0-p2), qualimap (v2.2) and preseq (v2.0.0). Differential gene expression was performed with DESeq2 package (v1.16.1, R v3.4.0) and with the same package read counts were normalised on the estimated size factors.

Resource       | URL
-------------- | --------------
GRCm38         | http://oct2016.archive.ensembl.org/index.html
FastQC         | http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
Trim_galore    | http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
TopHat2        | https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-4-r36
ClusterFlow    | http://dx.doi.org/10.12688/f1000research.10335.2
MultiQC        | http://dx.doi.org/10.1093/bioinformatics/btw354
HTSeq-counts   | http://dx.doi.org/10.1093/bioinformatics/btu638
Feature_counts | http://www.ncbi.nlm.nih.gov/pubmed/24227677
Qualimap       | https://doi.org/10.1093/bioinformatics/bts503
Preseq         | http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.2375.html

A custom module for TopHat2 double map is provided in this repository, and can be run, by copying it into the modules directory of a ClusterFlow installation. With HTSeq-Counts gene count tables figures in the table below can be reproduced with the R script provided in this repository.

### Script to reproduce paper figures ###

The provided R script assumes the script is placed in a directory containing a subdirectory (called HTSeq_Counts) with all the htseq-counts files (one per sample). The script can be run interactively in R-studio or as a batch using Rscript. Note that some of the figures in the manuscript have had some label positions moved manually to prevent overlaps.

Figure    | Description | Output Filename
--------- | ----------- | ------------------------
Figure 3A | PCA         | 2018-Balmas-Colucci_DESeq2_Annotated_PCA_Fig3A.pdf
Figure 3B | MA Plot     | 2018-Balmas-Colucci_DESeq2_MA_fc5_sig0.05_res_Uterine_Others_Fig3B.pdf
Figure 4A | Heatmap     | 2018-Balmas-Colucci_DESeq2_CountMatrixHeatmap_topDEGs_lf7.5sig0.05_Fig4A.pdf
Figure 4B | MA Plot     | 2018-Balmas-Colucci_DESeq2_MA_fc5_sig0.05_res_Virgin_E9.5_Fig4B.pdf
Figure 4D | MA Plot     | 2018-Balmas-Colucci_DESeq2_MA_fc5_sig0.05_res_Virgin_E18.5_Fig4D.pdf

Additional Data

Description                    | Output Filename
------------------------------ | ------------------------
DEG Results: Uterine_vs_Others | 2018-Balmas-Colucci_DESeq2_DEGs_Uterine_vs_Others.csv
DEG Results: Virgin_vs_E18.5   | 2018-Balmas-Colucci_DESeq2_DEGs_Virgin_vs_E18.5.csv
DEG Results: Virgin_vs_E9.5    | 2018-Balmas-Colucci_DESeq2_DEGs_Virgin_vs_E9.5.csv
Normalised reads counts        | 2018-Balmas-Colucci_DESeq2_NormalisedCounts.csv

### Sample Table ###


sampleName   | condition  | cell | tissue  
------------ | ---------- | ---- | -------  
SLX-10383.01 | Virgin     | ILC2 | Uterine  
SLX-10383.02 | Virgin     | ILC2 | Uterine  
SLX-10383.03 | Virgin     | ILC2 | Uterine  
SLX-10383.04 | E9.5       | ILC2 | Uterine  
SLX-10383.05 | E9.5       | ILC2 | Uterine  
SLX-10383.06 | E9.5       | ILC2 | Uterine  
SLX-10383.07 | E18.5      | ILC2 | Uterine  
SLX-10383.08 | E18.5      | ILC2 | Uterine  
SLX-10383.09 | E18.5      | ILC2 | Uterine  
SLX-9189.04  | Lymph_Node | ILC2 | Other    
SLX-9189.05  | Lung       | ILC2 | Other    
SLX-9189.06  | Lymph_Node | ILC2 | Other    
SLX-9189.07  | Lung       | ILC2 | Other    
SLX-9189.08  | Lymph_Node | ILC2 | Other    
SLX-9190.01  | Lung       | ILC2 | Other    


### Links ###

Description   | URL
------------- | ----------
Publication   | [bioRxiv]([bioRxiv](https://doi.org/10.1101/348755))
Raw Data      | ArrayExpress EMBL-EBI <br>Data to be released on publication<br>Uterine [E-MTAB-5803](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5803) <br> Lung/Lymph_Node [E-MTAB-5806](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5806)
Colucci Group | [Colucci group website](http://moffettcoluccilab.org/francesco-colucci/)

### Contact ###

Contact rsh46 -at- cam.ac.uk for bioinformatics related queries
