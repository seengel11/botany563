---
"Botany563 Project Notebook"
---

## Data 

Pheuticus is a genus of New World song birds also known as cardinals. Theses bird species are found across north, central, and south america. They are characterized by their complex songs, strong and large bills, and bright pigmentations. My data set contains information 22 Pheucticus samples and a single outgroup taxon. Sequenced are 9 loci (See table).

![](/Users/saraengel/Desktop/locus_table.png)

Data was taken from Pulgarin et. al:
https://www.sciencedirect.com/science/article/pii/S1055790313002297


The data was downloaded from the paperd Dryad repository:
https://datadryad.org/stash/dataset/doi:10.5061/dryad.844k8

### Quality Control
Aliments were conducted using Sequencher v4.9. Indelligent was used to resolve insertion/deletion events in introns in homologous nuclear alleles. They used the software PHASE v1.2 to infer the most likely arrangement of alleles at heterozygous sites in nuclear introns. Ambigious were labeled where a posterior probability was <0.90. Finally, They tested for recombination using RDP3. They constructed Gene Trees using MrBayes v3.1.2. The Akaike Information Criterion (AIC), as implemented in MrModeltest v2.2, was used to select the best model of nucleotide substitution for each locus.

## 2025-02-12
- Converted nexus files from paper to .fasta files for alignement using bugaco.com (https://sequenceconversion.bugaco.com/converter/biology/sequences/tab_to_fasta.php). Fiels downloaded to computer at Path:botany563/data/prealigned data. Files names follow this format 'locus'_seq.fasta.
- Aligned sequences using ClustalW (see 'MSA' for more info)

## MSA

### Clustal W

#### Description
Software:	
ClustalW

Description:
ClustalW is a widely-used tool for performing multiple sequence alignments of nucleic acids or proteins. It generates alignments in a hierarchical manner using a progressive method.
Strengths	
- User-friendly, widely accepted for multiple sequence alignment.
- Handles large datasets.
- Offers a range of output formats (e.g., Phylip, NEXUS).

Weaknesses	
- Can be slow with very large datasets.
- Doesn't always handle very divergent sequences well.
- Limited fine-tuning options compared to newer alignment tools.

Assumptions	
- Input sequences are assumed to be homologous.
- Assumes evolutionary relationships between the sequences.
- Assumes progressive alignment is the best approach.

User Choices	
- Choosing between different alignment modes (e.g., fast or accurate).
- Option to use various substitution matrices (e.g., PAM, BLOSUM).
- Setting gap penalties.

#### Terminal Commands: 
`cd saraengel/Desktop/botany563/data/Prealignment data`
`clustalw -ALIGN -INFILE='locus'_seq.fasta -OUTFILE='locus'-aligned.fasta -OUTPUT=FASTA`

#### File Location
- Output files located through following path: botany563/data/clustalw output
- Each locus alignment has to output files: 'locus'-aligned.fasta & 'locus'_seq.dnd

## Citation:

Pulgarín-R, P. C., Smith, B. T., Bryson, R. W., Spellman, G. M., & Klicka, J. (2013). Multilocus phylogeny and biogeography of the New World Pheucticus grosbeaks (Aves: Cardinalidae). Molecular Phylogenetics and Evolution, 69(3), 1222-1227. https://doi.org/10.1016/j.ympev.2013.05.022
