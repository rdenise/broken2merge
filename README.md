# Broken2merge

## Description
The `broken2merge` software is designed to concatenate the genes that belong to the same species but appear to be broken due to assembly issues. It provides a solution for merging fragmented gene sequences into a single, complete sequence.

## Features
- Gene concatenation: `broken2merge` identify and merge fragmented gene sequences into a single, continuous sequence.
- Assembly error detection: `broken2merge` includes error detection mechanisms to identify and handle assembly errors, ensuring accurate gene concatenation.

## Installation
To install `broken2merge`, follow these steps:

### Install with pipy:

```bash
pip install broken2merge
```

### Install in a conda/mamba env:

```bash
mamba create -n broken2merge python=3.12 tqdm biopython numpy
pip install broken2merge
```

## Usage

To use `broken2merge`:

```bash
usage: broken2merge [-h] [-v] [-V] -i FASTA_FILE [-o OUTPUT] [-s SEPARATOR]

Takes a alignement and find genes that seems broken and merge them together

options:
  -h, --help            show this help message and exit

General input dataset options:
  -v, --verbose         Show verbose output. (For debugging purposes)
  -V, --version         Show the version number and exit.
  -i FASTA_FILE, --input FASTA_FILE
                        Path to an input fasta file (Required)
  -o OUTPUT, --output OUTPUT
                        Path of the output folder (Default: merge_broken_res)
  -s SEPARATOR, --separator SEPARATOR
                        Separator to use to split the gene name (Default: ';')
  --force_merge         Force the merge of the genes even if they might be paralogs, will give a unaligned file as main input
```

## Example

Here's an example of how to use `broken2merge` to merge gene sequences:

```bash
broken2merge -i ftsK.aln.fas -o test -s ';'
```

Here for an input file named `ftsK.aln.fas` and output folder named `test` and the separator is `;` in the gene name (`species_name;gene_name`).


## Explaination of the concatenation

Case that will be handled:

1. Two genes that do not overlap in the alignment

```
speciesA_geneA  ------------------atgattgaactcgccc                
speciesA_geneA  atgattgaactcgccc------------------
```

it will become

```
speciesA_merge atgattgaactcgcc--catgattgaactcgccc
```

2. Two genes that overlap and perfectly in the alignment in only one part

```
speciesA_geneA  ----------------actcgcccatgattgaactcgccc
                                ||||||||
speciesA_geneB  atgattgaactcgcccactcgccc----------------
```

will become

```
speciesA_merge atgattgaactcgcccactcgcccatgattgaactcgccc
```

3. Two genes with a overlap not at the extremity of the gene

```
speciesA_geneA  ---------------------actcgcccatgattgaactcgccc
                                     ||||||||
speciesA_geneB  atgattgaactcgccc-----actcgccc----------------
```


will become

```
speciesA_merge atgattgaactcgccc-----actcgcccatgattgaactcgccc
```

4. Two genes wit an overlap that is not 100% perfect

```
speciesA_geneA  ---------------------actcgcccatgattgaactcgccc
                                     || ||| | 
speciesA_geneB  atgattgaactcgccc-----acccgcgc---------------
```

The sequences will be discarded from the alignement.

If `--force_merge` is used, it the sequences will be merged with the sequence one after the other in order from the alignement.
The output fasta file will unaligned.

it will become

```
speciesA_merge atgattgaactcgcccacccgcgcactcgcccatgattgaactcgccc
```


