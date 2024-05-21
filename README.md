# rnaCrosslinkOO


The rnaCrosslinkOO R package is designed to analyse data from RNA crosslinking and high-throughput sequencing experiments. The package has 3 main steps, clustering, cluster trimming and folding. Raw and processed data is stored at each stage in the analysis and plot and metrics can be obtained through the use of the packages other functions. 


A vignette is available at: [vignette](https://cran.r-project.org/package=rnaCrosslinkOO). It is recommended to install the package through CRAN. 

Datasets can be found in the [paper](https://academic.oup.com/bioinformatics/article/40/4/btae193/7643507) but we often have unpublished in house datasets that we are willing to share. So please get in touch. 

```
# CRANs
install.packages("rnaCrosslinkOO")
# GITHUB
devtools::install_github("JLP-BioInf/rnaCrosslinkOO")
```


# Input for rnaCrosslinkOO
\
The main input files for rnaCrosslinkOO is a tab delimited text file containing the reads and mapping location on the transcriptome. This can be manually created if your library preparation protocol does not suit the pre-processing pipeline. Although, the easiest way to obtain these files is to use the Nextflow pipeline   [here](https://github.com/JLP-BioInf/rnaCrosslinkNF), there are also details in the  [vignette](https://cran.r-project.org/package=rnaCrosslinkOO). There is test data that ships with the package, this contains data for the 18S rRNA and it's interactions with the 28S rRNA. However, full data-set already published can be found here:[Un-enriched rRNA dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246412), the [SARS-CoV-2 dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154662).

Pre-requisites:

1. Install the rnaCrosslinkOO package
2. Input files (nexflow, custom or downloaded)
3. Meta-data table
4. ID of the RNA of interest (from the transcript reference )
5. A fasta sequence of the RNA of interest  (from the transcript reference )
6. Locally install [The vienna package](https://www.tbi.univie.ac.at/RNA/)
7. A set of interactions to compare to (optional)
8. Reactivities (optional)


# The COMRADES experiment
\
The COMRADES (crosslinking of Matched RNA and Deep Sequencing) experimental protocol for the prediction of RNA structure in vivo was first published in 2018 (Ziv et al., 2019). The protocol has subsequently been use to analyse the structure of SARS-CoV-2 (Ziv et al., 2020) and NORAD:

* COMRADES determines in vivo RNA structures and interactions. (2018). Omer Ziv, Marta Gabryelska, Aaron Lun, Luca Gebert. Jessica Sheu-Gruttadauria and Luke Meredith, Zhong-Yu Liu,  Chun Kit Kwok, Cheng-Feng Qin, Ian MacRae, Ian Goodfellow , John Marioni, Grzegorz Kudla, Eric Miska.  Nature Methods. Volume 15. https://doi.org/10.1038/s41592-018-0121-0   

* The Short- and Long-Range RNA-RNA Interactome of SARS-CoV-2. (2020). Omer Ziv, Jonathan Price, Lyudmila Shalamova, Tsveta Kamenova, Ian Goodfellow, Friedemann Weber, Eric A. Miska. Molecular Cell,
Volume 80. https://doi.org/10.1016/j.molcel.2020.11.004


Ziv et al., 2020. "Virus-inoculated cells are crosslinked using clickable psoralen. Viral RNA is pulled down from the cell lysate using an array of biotinylated DNA probes, following digestion of the DNA probes and fragmentation of the RNA. Biotin is attached to crosslinked RNA duplexes via click chemistry, enabling pulling down crosslinked RNA using streptavidin beads. Half of the RNA duplexes are proximity-ligated, following reversal of the crosslinking to enable sequencing. The other half serves as a control, in which crosslink reversal proceeds the proximity ligation. "

After sequencing, short reads are produced similar to a spliced / chimeric RNA read but where one half of the read corresponds to one half of a structural RNA duplex and the other half of the reads corresponds to the other half of the structural RNA duplex. This package has been designed to analyse this data. Although the package was designed around the analysis of COMRADES data, the package will accept data provided it is the right format. 

---
