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


## 3/1/25
- Created guide trees using distance and parsimony methods of A40-aligned.fasta
- Trees created in a .Rmd document located in /Users/saraengel/Desktop/botany563/scripts labeld Trees.Rmd

## Distance and Parsimony
### R Packages:
```{r}
library(ape)
library(adegenet)
library(phangorn)
```

### Distance Based Methods

 Neighbour Joining

Description:
The Neighbor-Joining (NJ) method is a distance-based algorithm used to construct phylogenetic trees. It works by clustering sequences based on pairwise genetic distances and aims to minimize the total branch length of the resulting tree. The NJ method is popular for its simplicity and speed, especially when handling large datasets.

Strengths:
- Efficiency: The NJ algorithm is computationally efficient, making it suitable for large datasets.
- Speed: It is faster than many other tree construction methods (e.g., Maximum Likelihood or Bayesian inference) since it only requires distance calculations between sequences.
- Simplicity: The NJ method is easy to understand and implement, requiring only a distance matrix as input.
- Applicability: It works well for cases where you have reliable distance measures between taxa and is often used in molecular biology and evolutionary studies.

Weaknesses:
- Assumes Additivity: The NJ method assumes that the evolutionary distances between species are additive (i.e., the total distance between two taxa is the sum of the distances through the internal nodes), which may not hold true in all cases.
- Loss of Information: Distance-based methods like NJ can lose some information compared to character-based methods (e.g., Maximum Likelihood), as they reduce sequence data into pairwise distances.
- Susceptible to Errors in Distance Estimation: The accuracy of the tree depends on the quality of the distance matrix. Any errors in estimating genetic distances can propagate to errors in the tree structure.
 - Lacks Statistical Support: The NJ method does not provide a statistical measure of confidence (e.g., posterior probabilities or bootstrap values) for the branches, which are often important in assessing tree reliability.

Assumptions:
- Distance Measure: The accuracy of the NJ tree heavily depends on the assumption that the genetic distance measure used (e.g., Jukes-Cantor, Kimura, etc.) accurately reflects evolutionary relationships. 
- Additivity: The method assumes that evolutionary distances are additive, which may not always be valid, particularly in cases with horizontal gene transfer or varying rates of evolution across lineages.

####R Script
```{r}
dna <- fasta2DNAbin(file="~/Desktop/botany563/data/clustalw output/A40-aligned.fasta")
D <- dist.dna(dna, model="TN93")
tre <- nj(D)
tre <- ladderize(tre)
plot(tre, cex=.6)
title("A simple NJ tree")
```
### Parsimony Methods
Tree building using Parsmimony: 

Description:
The parsimony method is a character-based approach for constructing phylogenetic trees. It selects the tree topology that minimizes the total number of evolutionary changes (e.g., mutations, insertions, deletions) required to explain the sequence data. Parsimony is based on the principle that the simplest explanation, or the one with the least change, is the most likely.

Strengths:
- Simplicity: Parsimony methods are easy to understand and interpret. They directly reflect the principle of Occam's Razor—preferring the simplest solution that requires the fewest evolutionary events.
- Computational Feasibility: Parsimony can be computationally efficient for smaller datasets or fewer taxa.
- haracter-Based: Unlike distance-based methods, parsimony considers individual characters (nucleotide or amino acid positions), retaining more information from the original sequence data.
- No Assumption of Molecular Clock: Parsimony does not require assuming a constant rate of evolution across lineages, which can make it more flexible for certain types of data.

Weaknesses
- Inconsistent with Long Branch Attraction: Parsimony methods can be prone to "long branch attraction," where distantly related taxa with high rates of evolution appear artificially close due to many convergent changes.
- Computational Complexity with Large Datasets: As the number of taxa and sequences grows, the number of possible tree topologies increases exponentially, making parsimony searches computationally intensive.
- No Model of Evolution: Parsimony lacks an explicit model of sequence evolution. This can be a weakness because it does not account for unequal rates of mutation at different sites or among different lineages.
- Limited Statistical Support: Like NJ, the parsimony method does not inherently provide statistical support (e.g., posterior probabilities) for the branches, which is often crucial for assessing tree reliability.

Assumptions
- Character Independence: Parsimony assumes that each character (nucleotide or amino acid site) evolves independently of others.
- Equal Weighting of Changes: Parsimony assumes that all character changes are equally probable. Some variants of the method allow for differential weighting of transitions, transversions, or indels, but standard parsimony assigns equal weights.
- Minimal Evolution: Parsimony assumes that the simplest tree (requiring the fewest changes) is the most likely explanation for the observed data, without considering varying rates of evolution across lineages.

#### R Script

```{r}
dna <- fasta2DNAbin(file="~/Desktop/botany563/data/clustalw output/A40-aligned.fasta")
dna2 <- as.phyDat(dna)
tre.ini <- nj(dist.dna(dna,model="raw"))
parsimony(tre.ini, dna2)
tre.pars <- optim.parsimony(tre.ini, dna2)

plot(tre.pars, cex=0.6)
```

## 3/19/25
- Used IQ-tree to make maximum likelihood tree of A40 gene. 

### Maximum Likelihood

IQ Tree
Software Description:	
IQ-TREE is a fast and efficient software for inferring phylogenetic trees using maximum likelihood (ML) methods. It incorporates model selection, ultrafast bootstrapping, and other advanced features to ensure accurate and robust tree estimation. It's widely used in molecular phylogenetics, especially for large datasets.

Strengths	
- Fast and scalable: Optimized for large datasets with high performance.
- Ultrafast bootstrap: Provides rapid, accurate bootstrap support values (UFBoot2).
- Automated model selection: Automatically selects the best-fit substitution model.
- Bayesian-like branch supports: SH-like approximate likelihood ratio test (SH-aLRT) for branch support.
- Wide range of models: Supports nucleotide, amino acid, and codon models for sequence evolution.

Weaknesses
- ML focus: IQ-TREE is focused on maximum likelihood methods and may not be the best choice if you're looking for Bayesian methods or non-parametric techniques.
- Complexity: Advanced features require familiarity with phylogenetic methods, making it harder for beginners.
- Assumes independence: As with most phylogenetic models, it assumes independence between sites, which may not hold for all datasets.

Assumptions	
- Evolutionary model assumptions: Assumes a substitution model for sequence evolution (e.g., GTR, HKY85, etc.), which may not always perfectly fit the data.
- Independence of sites: Assumes that nucleotide or amino acid sites evolve independently.
- Tree-like evolution: Assumes the data can be represented by a bifurcating tree.
- Stationarity: Assumes that the substitution process is stationary (i.e., unchanged across the tree).

User Choices	
- Substitution model: Users can specify a substitution model or allow IQ-TREE to select the best model.
- Bootstrap methods: Choose between ultrafast bootstrap or other support estimation methods.
- Partitioning schemes: Ability to partition data and apply different models to different parts of the dataset.
- Branch support metrics: Option to use SH-aLRT, ultrafast bootstrap, or other branch support values.
- Gene-wise partitioning: Can handle gene partitions or whole-genome alignments with model selection for each partition.


#### Script

##### Terminal for ML
 
```
cd Desktop/iqtree-1.6.12-MacOSX
bin\iqtree -s /Users/saraengel/Desktop/botany563/data/clustalw output/A40_seq.dnd

```
File outputs where placed in the following path: /Users/saraengel/Desktop/botany563/data/Maximum Likelihood Output/A40, ML
 
##### R For plotting
Rooted tree at the outgroup(Granatellus venusta, GR JK04-078) as stated by the original paper (Pulgarin et al.)
Text was copied from A40-aligned.fasta.iqtree file from IQ-tree output
```{r}

tree <- read.tree(text="(Pms_BRB1162a:0.0000020199,Pmw_DHB5399a:0.0000028651,((Pme_dhb4768b:0.0024695124,((((((Pc_116284a:0.0000022538,Pc_116283a:0.0000020689):0.0019400077,(Pac_TA22a:0.0000020646,Pac_TA22b:0.0000025308):0.0000020207):0.0000020423,((Pa_B8248a:0.0000028651,Pa_B8248b:0.0000020667):0.0000020443,Pc_B5186b:0.0019399415):0.0019399892):0.0040621032,(((Py_BTS08324a:0.0000020249,Py_BTS08324b:0.0000020292):0.0019393181,Ps_SIN152b:0.0000025567):0.0019396388,((Ps_SIN152a:0.0000020423,Py_BTS08305a:0.0000022396):0.0000020596,(Py_BTS08305b:0.0039327275,(Pmw_DHB5399b:0.0000024018,Pme_JK09429b:0.0000020575):0.0039376824):0.0000022921):0.0000020247):0.0040694295):0.0019357374,(Gr_JK04078a:0.0000024023,Gr_JK04078b:0.0019393843):0.0129577583):0.0000020075,(((((Pl_DHB5709a:0.0000028651,Pl_DHB5709b:0.0000028651):0.0000029230,Pl_TKA28b:0.0000028651):0.0000022724,Pl_TKA28a:0.0000028651):0.0000021622,Pl_DHB5227a:0.0000020429):0.0019695850,((Pt_B16050a:0.0000020966,Pt_B16050b:0.0000021003):0.0041204754,Pl_DHB5227b:0.0019617880):0.0020001200):0.0019739723):0.0019814287):0.0000029503,((Pms_BRB1162b:0.0019757602,Pme_JK09429a:0.0019733915):0.0017116462,Pme_dhb4768a:0.0082756970):0.0065230974):0.0000020831);")


plot(tree)
nodelabels()

rooted_tree = root(tree, node = 51, resolve.root = TRUE)
plot(rooted_tree)
```


## Citation:

Pulgarín-R, P. C., Smith, B. T., Bryson, R. W., Spellman, G. M., & Klicka, J. (2013). Multilocus phylogeny and biogeography of the New World Pheucticus grosbeaks (Aves: Cardinalidae). Molecular Phylogenetics and Evolution, 69(3), 1222-1227. https://doi.org/10.1016/j.ympev.2013.05.022
