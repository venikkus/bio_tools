# Bio tools

**Bio tools** is a set of packages for filtering FASTQ files, processing DNA or RNA sequences and some types of bioinformatics data like ``.gbk``, 

---
Authors:
* **Software:** [venikkus](https://github.com/venikkus), Saint-Petersburg, Russia.

* **Idea, supervisor:** [Bioinformatics Institute](https://bioinf.me/en), Saint-Petersburg, Russia.

---
## Content

- [Bio tools](#bio-tools)
- [Content](#content)
- [Install](#install)
- [Running instructions](#running-instructions)
- [Features](#features)
- [Examples](#examples)
- [Contact](#contact)

## Install

You need to download the folder to your PC. You can download directly or
clone the repository to your computer.

```bash
git clone git@github.com:venikkus/bio_tools.git
```
---
## Running instructions

To work with the package, you can import the main script (bio_tools or bio_files_processor) and call any of the functions.

## Features

There are two scripts, which have 5 functions and additional modules for them.
### 1. bio_tools

This script has 2 main function and additional additional modules for them.
   
   1. **run_dna_rna_tools** takes a string and an operation and performs one of actions.

   - transcribe — transcribes a sequence
   - complement — converts a forward sequence to a complementary
   - reverse — reverses a sequence
   - reverse_complement — converts a forward sequence to a reverse complementary sequence
   - annealing_temperature — calculates the annealing temperature of a sequence
   - is_palindrome — checks the supplied sequence is a biopalindrome
   - is_nucleotide — checks the supplied sequence is a nucleotide
   - is_primer — checks the supplied sequence is a primer

   2. **filter_fastq** takes raw ``.fasta`` file with reads and reading quality indicators and calculates GC composition of the read, the quality of the read and save filtered data in ``.txt`` file. The module contains auxiliary functions that help with output, file recording, and calculation of filter parameters:
   - calculate_gc_bounds — calculates the GC content of a nucleotide sequence
   - calculate_quality_threshold — calculates the average Q-score of a nucleotide sequence from a FASTQ file
   - make_bounds — creates bounds from an integer or returns the given bounds as is
   - is_bounded — checks if a value is within the specified bounds
   - read_fastq — reads a FASTQ file and converts it to a dictionary
   - write_fastq — writes the filtered sequences to a FASTQ file




### 2. bio_files_processor

   1. **convert_multiline_fasta_to_oneline** — converts a multi-line FASTA file to a single-line FASTA format.

   This function reads the FASTA file where sequences might be spread across
   multiple lines and outputs the new FASTA file or prints the sequences, where
   each sequence is written with single line.


   2. **select_genes_from_gbk_to_fasta** — parses the output of the BLAST search and extracts the first description line
    for each query.

   This function reads the BLAST results file, extracts descriptions
   for each query, and outputs sorted results to file or prints
   them to console.

   3. **select_genes_from_gbk_to_fasta** extracts genes around gene of interest.

   Extracts specified gene sequences along with n_before upstream
   and n_after downstream genes from a GenBank file.

All script functions collect data into the file if you specify path to it. Otherwise, result is output to console.

## Examples

Bio_tools contains 8 simple function. For example, let's try transcribe function:
```
run_dna_rna_tools("AGt", "transcribe")
```
Result of the function is transcribed sequence:
```
"AGu"
```

Let`s try **bio_files_processor** functions with our [example data](https://github.com/venikkus/bio_tools/tree/add_bio_tools/data).


```
convert_multiline_fasta_to_oneline(
    input_fasta="data/example_multiline_fasta.fasta",
    output_fasta="data/oneline_fasta.fasta",
)
```

As we can see in the output, every sequence is concatenate to the one line (only 3 examples shown):
```
>5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+)
ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCTGGTGCTGTG
>16S_rRNA::NODE_80_length_720_cov_1.094737:313-719(+)
TTGGCTTCTTAGAGGGACTTTTGATGTTTAATCAAAGGAAGTTTGAGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGACAAAGTCAACGAGTTTTATTATTATTCCTTTATTGAAAAATATGGGTAATCTTGTTAAACTTTGTCGTGCTGGGGATAGAGCATTGCAATTATTGCTCTTCAACGAGGAATTCCTAGTAAGCGCAAGTCATCAGCTTGCGTTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGAATGGCTTAGTGAGCCCTTGGGAGTGGTCCATTTGAGCCGGCAACGGCACGTTTGGACTGCAAACTTGGGCAAACTTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGT
>16S_rRNA::NODE_1_length_2558431_cov_75.185164:2153860-2155398(+)
TTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTTGTTGGTGGGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTT
...
```

The function accepts the result of the BLAST algorithm.

```
parse_blast_output(
    input_file="data/example_blast_results.txt",
    output_file="data/blast_output.txt",
)
```

``blast_output.txt`` file contains 37 string with gene names from Description column in input file in alphabetical order (only 10 examples shown):
```
DNA methylase [Enterobacteriaceae]
DUF1380 domain-containing protein [Escherichia coli]
DUF1380 family protein [Enterobacteriaceae]
DUF4158 domain-containing protein [Klebsiella pneumoniae]
DUF905 domain-containing protein [Shigella sonnei]
DinI-like family protein [Escherichia coli]
KlcA [Escherichia coli]
PilI type IV pilus biogenesis protein [Salmonella enterica]
PilK [Escherichia coli]
PilN family type IVB pilus formation outer membrane protein...
...
```

File ``.gbk`` stores DNA and protein sequences. This function searches for genes of interest and those located near them.

```
select_genes_from_gbk_to_fasta(
    "data/example_gbk.gbk",
    ["ybcO"],
    n_before=2,
    n_after=2,
    output_fasta="data/TEST2_output_fasta.fasta",
)
```

As a result of the function's operation, ``.fasta`` file is obtained with the following content (only 1 example shown):
```
xerC_1
MGRRRSHERRDLPPNLYIRNNGYYCYRDPRTGKEFGLGRDRRIAITEAIQANIELFSGHKHKPLTARINSDNSVTLHSWLDRYEKILASRGIKQKTLINYMSKIKAIRRGLPDAPLEDITTKEIAAMLNGYIDEGKAASAKLIRSTLSDAFREAIAEGHITTNPVAATRAAKSEVRRSRLTADEYLKIYQAAESSPCWLRLAMELAVVTGQRVGDLCEMKWSDIVDGYLYVEQSKTGVKIAIPTALHVDALGISMKETLDKCKEILGGETIIASTRREPLSSGTVSRYFMRARKASGLSFEGDPPTFHELRSLSARLYEKQISDKFAQHLLGHK

>emrE_1
MNPYIYLGGAILAEVIGTTLMKFSEGFTRLWPSVGTIICYCASF

>ybcO_1 # gene of interest
MADLRKAARSRECQVRIPGVCNGNPETSVLAHIRLTGLCGTGTK

>rusA_1
MLDIGLAMPVKIRIECHMPDRRRRDLDNLQKAAFDALTKAGFWL

>ylcG
MFEFYMAERLRHRWGRLRLYRFPGSVLTDYRILKNYAKTLTGAG
...
```



## Contact

Please report any problems directly to the GitHub
[issue tracker](https://github.com/venikkus/bio_tools/issues).<br/>
Also, you can send your feedback to
[niksamusik@gmail.com](mailto:niksamusik@gmail.com).
