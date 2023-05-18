# Homology Evaluation of Xenopeptides (HEX)

## Background

HEX is designed to identify tumor peptides with high similarity to viral or bacterial epitopes in order to find candidate peptides that perform molecular mimicry. This has been explored in autoimmunity and HEX can extend that into anti-tumour responses.

## Installation
To run the installation, open the ```server.R``` file in Rstudio and press "run app" in the top right hand corner. Rstudio should have a prompt to install any of the necessary packages that aren't installed yet. Alterntively, visit https://whalleyt.shinyapps.io/homology-evalutation-xenopeptides/ to run the tool online.

## Usage
There are two boxes for inputting sequences, one for the reference sequence (typically tumour proteins) and query peptides (typically foreign protein sequences coming from a bacterial or viral species of interest).
The user must upload these sequences in FASTA format (see below). The sequences can be uploaded either as a file or sequence. When running locally, one can use the pre-compiled database however this is disabled on the web app due to resource limitations.

The user then can select their method from MHC affinity prediction (see [here](http://tools.iedb.org/mhci/help/) for references and discussion about the different methods), epitope length of interest, aligment matrix and organism/MHC allele.

An example FASTA file format is below:

```
>sequence 1
THISISASEQUENCE
>sequence 2
THISISANOTHERSEQ

```

The output is displayed in a table which can filtered. The ```get best``` option sorts them by IC50, ```by alignment``` by sequence homology and ```by ALS score``` by the pseudo-combinatorial peptide library score. *Please note, if you are working with a low number of peptides it is sensible to remove the IC50 and alignment score threshold as no peptides may satisfy this filter.

## Citation
[Chiaro J, Kasanen HHE, Whalley T, Capasso C, Grönholm M, Feola S, Peltonen K, Hamdan F, Hernberg M, Mäkelä S, Karhapää H, Brown PE, Martins B, Fusciello M, Ylösmäki EO, Greco D, Kreutzman AS, Mustjoki S, Szomolay B, Cerullo V. Viral Molecular Mimicry Influences the Antitumor Immune Response in Murine and Human Melanoma. Cancer Immunol Res. 2021 Aug;9(8):981-993. doi: 10.1158/2326-6066.CIR-20-0814. Epub 2021 Jun 8. PMID: 34103348; PMCID: PMC8974425.](https://aacrjournals.org/cancerres/article/81/13_Supplement/1488/667320/Abstract-1488-Viral-molecular-mimicry-influences)
