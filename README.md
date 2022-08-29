# scDoRI
Save progress of my labrotation

The analysis of this project uses a multiome mouse gastrualation dataset from (Argelaguet, bioRxiv, 2022) doi: https://doi.org/10.1101/2022.06.15.496239.


* `ArchK` contains R scripts with adapted ArchR (www.ArchRProject.com) functions. Cloning the ArchR repository (https://github.com/GreenleafLab/ArchR) and adding these R scripts as well as the function names to the name space enables the user to use the new functions.
* `Rmds` contains Rmarkdowns with the early analysis of the mouse gastrulation dataset, including 4 timepoints, but no perturbation.
* `Rmds_perturbation contains the same analysis as above, but for the complete dataset, including knockout of the transcription factor brachyury.
* `job_scripts` contains bash scripts as well as R and python scripts that were run on the Cluster, because they require a lot of memory or a GPU.
* `jupyter` contains any analysis that was conducted in Python. 
* `old` contains outdated analysis scripts which might still be informative.


**Abstract of the project**

Understanding the individual players within gene regulatory networks, meaning which transcription factors cause transcription or repression of a gene, is of crucial importance for understanding how different cell types and tissues are established, but also how abnormal gene regulation due to mutations can impair development and cause disease. Single-cell multiome data in combination with computational tools developed in recent years, enable researchers to make predictions about gene regulatory links. It remains difficult to validate predictions of these tools, because biological experiments are time-consuming and expensive. Here, we compute and compare gene activity scores based on peak-to-gene links inferred with ArchR and try to validate them based on the correlation between gene activity scores and gene expression. This approach aids the understanding and interpretation of gene activity scores and regulatory links. A known shortcoming of transcription factor binding motifs from public databases is that transcription factors of the same family often have the same motif, even though \textit{in vivo} they bind to different regions of the genome. We show that transcription factor binding probabilities predicted using convolutional neural networks do not improve the distinction between ChromVar motif enrichment scores of transcription factors of the same family, meaning new tools are required to improve this distinction. Pseudotime trajectories  were computed for the Erythroid and Endothelium lineage to demonstrate a temporal delay between chromatin accessibility and gene expression. Pseudotime for other cell lineages was not in concordance with biological knowledge, highlighting the limitations of these tools and the need to improve trajectory models further. This work highlights some shortcomings of available tools, as well as, the importance of biological validation of computational predictions. 
