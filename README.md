# comradesOO 

# COMRADES experiment


The COMRADES experimental protocol for the prediction of RNA structure in vivo was first published in 2018 (Ziv et al., 2019) where they predicted the structure of the Zika virus. The protocol has subsequently been use to predict the structure of SARS-CoV-2 (Ziv et al., 2020). Have a look to get an understanding of the protocol:

* COMRADES determines in vivo RNA structures and interactions. (2018). Omer Ziv, Marta Gabryelska, Aaron Lun, Luca Gebert. Jessica Sheu-Gruttadauria and Luke Meredith, Zhong-Yu Liu,  Chun Kit Kwok, Cheng-Feng Qin, Ian MacRae, Ian Goodfellow , John Marioni, Grzegorz Kudla, Eric Miska.  Nature Methods. Volume 15. https://doi.org/10.1038/s41592-018-0121-0   

* The Short- and Long-Range RNA-RNA Interactome of SARS-CoV-2. (2020). Omer Ziv, Jonathan Price, Lyudmila Shalamova, Tsveta Kamenova, Ian Goodfellow, Friedemann Weber, Eric A. Miska. Molecular Cell,
Volume 80
    https://doi.org/10.1016/j.molcel.2020.11.004
    


![Figure from Ziv et al., 2020. Virus-inoculated cells are crosslinked using clickable psoralen. Viral RNA is pulled down from the cell lysate using an array of biotinylated DNA probes, following digestion of the DNA probes and fragmentation of the RNA. Biotin is attached to crosslinked RNA duplexes via click chemistry, enabling pulling down crosslinked RNA using streptavidin beads. Half of the RNA duplexes are proximity-ligated, following reversal of the crosslinking to enable sequencing. The other half serves as a control, in which crosslink reversal proceeds the proximity ligation](https://github.com/JLP-BioInf/comradesOO/tree/master/vignettes/comradesProtocol.jpg)


After sequencing, short reads are produced where one half of the read corresponds to one half of an RNA duplex and the other half of the reads corresponds to the other half of the RNA duplex. This package has been designed to analyse this data. The short reads need to be processed is a specific way, see the next section. 




---

# COMRADES data pre-processing

## Nextflow pipeline 
\
Fastq files produced from the comrades experiment can be processed for input into the comradesOO using the Nextflow pre-processing pipeline, to get more information visit here. https://github.com/JLP-BioInf/comradesnf . The pipeline Takes the reads through trimming alignment, QC and the production of the files necessary for input to comradesOO.



## Nextflow pipeline output
\
The main output files are the files entitled *X_gapped.txt*. These are the input files for comradesOO. The columns of the output files are as follows:

You can create your own versions of this file if you have data using different library protocol. As long as you have 1 file for each sample with the following columns:

1. Read Name
2. Read Sequence
3. Side 1 transcript ID
4. Side 1 Position start in read sequence
5. Side 1 Position end in read sequence
6. Side 1 Coordinate start in transcript
7. Side 1 Coordinate end in transcript
8. NA
9. Side 1 transcript ID
10. Side 1 Position start in read sequence
11. Side 1 Position end in read sequence
12. Side 1 Coordinate start in transcript
13. Side 1 Coordinate end in transcript
14. NA





---

# Input for comrades-OO
\
The main input files for comrades-OO is a tab delimited text file containing the reads and mapping location on the transcriptome. This can be manually created although the easiest way to obtain these files is to use the nextflow pipeline detailed above. There is test data that ships with the package, this contains data for the 18S rRNA and it's interactions with the 28S rRNA. However, full data-sets already published can be found here:  [SARS-CoV-2 Dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154662) (files ending in "txt.gz") 


Pre-requisites:

1. Install the comradesOO package
2. Input files (nexflow, custom or downloaded)
3. Meta-data table
4. ID of the RNA of interest (from the transcript reference )
5. A fasta sequence of the RNA of interest  (from the transcript reference )
6. A set of interactions to compare to (optional)
6. Reactivities (optional)

See the package Vignette for full instructions 




