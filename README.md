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

cd bio_tools
pip install -r requirements.txt
```
---
## Running instructions

To work with the package, you can import the main script (bio_tools or bio_files_processor) and call any of the functions.

## Features

### 1. **bio_tools.py**
This script provides core functionality for processing nucleic acid and protein sequences:

#### **Biological Sequences Processing**
- **NucleicAcidSequence** — class for RNA and DNA sequences.
  - `.get_complement()` - complementing the sequence.
  - `.reverse()` - reversing the sequence.
  - `.reverse_complement()` - obtaining the reverse complement.
  - `.annealing_temperature()`- Calculating the annealing temperature.
  - `.check_palindrome()`- Checking if the sequence is a valid primer or palindrome.

- **DNASequence** — class for DNA sequences.
  - `.transcribe()` — transcribes DNA to RNA.
  - `.get_complement()` — computes the complementary sequence.
  - `.reverse_complement()` — computes the reverse complement.
  - `.annealing_temperature()` — calculates the annealing temperature.
  - `.check_primer()` — determines if the sequence is a valid primer.
  - `.check_palindrome()` — checks if the sequence is a palindrome.

- **RNASequence** — class for RNA sequences.
  - `.reverse_transcribe()` — Converts RNA back into DNA.
  - `.reverse_complement()` — Computes the reverse complement.

- **AminoAcidSequence** — class for protein sequences.
  - `.molecular_weight()` — Calculates the molecular weight of a protein.

#### **FASTQ Processing** 
- **filter_fastq()** — takes raw ``.fasta`` file with reads and reading quality indicators and calculates GC composition of the read, the quality of the read and save filtered data in ``.txt`` file. Filters a FASTQ file based on:
  - GC content,
  - Sequence length,
  - Quality score.

### 2. bio_files_processor
`Some addictionaly functions.`

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

Class DNASequence:
```python
dna = DNASequence("ATGCGT")
print(f"DNA: {dna}")
print(f"Complement: {dna.get_complement()}")
print(f"Reverse complement: {dna.reverse_complement()}")
print(f"Transcribed RNA: {dna.transcribe()}")
```

Output:
```python
DNA: ATGCGT
Complement: TACGCA
Reverse complement: ACGCAT
Transcribed RNA: AUGCGU
```

Class RNASequence:
```python
rna = RNASequence("AUGCGU")
print(f"RNA: {rna}")
print(f"Complement: {rna.get_complement()}")
print(f"Reverse complement: {rna.reverse_complement()}")
print(f"Reverse transcribed DNA: {rna.reverse_transcribe()}")
```
Output:
```python
RNA: AUGCGU
Complement: UACGCA
Reverse complement: ACGCAU
Reverse transcribed DNA: ATGCGT
```

Class AminoAcidSequence:
```python
protein = AminoAcidSequence("MVK")
print(f"Protein: {protein}")
print(f"Molecular weight: {protein.molecular_weight()}")
```

Output:
```python
Protein: MVK
Molecular weight: 412.5
```

If you try on DNASequence or RNASequence on invalid symbols or blank line into sequences it will raise custom error:

```python
DNASequence("ATGX")
DNASequence("")

FAILED dna_rna_tools_test.py::test_invalid_dna - ValueError: Invalid symbols: {'X'}
FAILED dna_rna_tools_test.py::test_empty_string - ValueError: Sequence cannot be empty.
```

To use `filter_fastq` function you need define input/output paths and set GC and length bounds, quality threshold.
```python
filter_fastq(
    input_fastq="data/example_fastq.fastq",
    output_fastq="filtered/filtered_output.fastq",
    gc_bounds=(40, 60),
    length_bounds=(20, 35),
    quality_threshold=30
)
```
Output:
```python
@SRX079804:1:SRR292678:1:1101:829239:829239 1:N:0:1 BH:ok
TCGATCCTTCTGCCTCAAAGTATACTAGGACGCAT
+
GGGDFGGBGFFEBFEDCBCDCGGGGBEEE=GE?EE
@SRX079804:1:SRR292678:1:1101:868419:868419 1:N:0:1 BH:ok
ATTCGTCAGGCCCAATAACATCATGAATTTCCAG
+
DEEEEEEEBDFFFFFFFF8FEED8@FFFBFFEFF
@SRX079804:1:SRR292678:1:1101:918742:918742 1:N:0:1 BH:failed
CTCTCCATGCACAAAGAATATCACAGCCAAA
+
EEEBA?@;B@EEE@BEE=?EDDDDADCDA?E
@SRX079804:1:SRR292678:1:1101:933189:933189 1:N:0:1 BH:failed
GTCTGCACTATCGAGGGCTGTGCCTTTGC
+
FEFFDBFF8FE>?DFFFCEBCEEBBEDE6
@SRX079804:1:SRR292678:1:1101:940351:940351 1:N:0:1 BH:changed:1
TGCCGTGGGAATGACAAACAAGCATCC
+
DECC@GFFBF=EBEAFDFGD?FFF8FF
@SRX079804:1:SRR292678:1:1101:940693:940693 1:N:0:1 BH:failed
CACATTATGAACTATGGGCACTGCAT
+
EEEGFDEDFEGGGGGFEGBGGGFGGG
@SRX079804:1:SRR292678:1:1101:955819:955819 1:N:0:1 BH:failed
CACCTAGCAGCAACGGACGAGTCAG
+
GGGGGEEEGGEGGGFGEGG;F@EFF
@SRX079804:1:SRR292678:1:1101:996098:996098 1:N:0:1 BH:failed
CTAAGAGAGTTTGTAATGCGGAC
+
DD=DBDBDC4EFFFD@?CD@ACD
```


Let`s try **bio_files_processor** functions with our [example data](https://github.com/venikkus/bio_tools/tree/add_bio_tools/data).


```python
convert_multiline_fasta_to_oneline(
    input_fasta="data/example_multiline_fasta.fasta",
    output_fasta="data/oneline_fasta.fasta",
)
```

As we can see in the output, every sequence is concatenate to the one line (only 3 examples shown):
```python
>5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+)
ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCTGGTGCTGTG
>16S_rRNA::NODE_80_length_720_cov_1.094737:313-719(+)
TTGGCTTCTTAGAGGGACTTTTGATGTTTAATCAAAGGAAGTTTGAGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGACAAAGTCAACGAGTTTTATTATTATTCCTTTATTGAAAAATATGGGTAATCTTGTTAAACTTTGTCGTGCTGGGGATAGAGCATTGCAATTATTGCTCTTCAACGAGGAATTCCTAGTAAGCGCAAGTCATCAGCTTGCGTTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGAATGGCTTAGTGAGCCCTTGGGAGTGGTCCATTTGAGCCGGCAACGGCACGTTTGGACTGCAAACTTGGGCAAACTTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGT
>16S_rRNA::NODE_1_length_2558431_cov_75.185164:2153860-2155398(+)
TTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTTGTTGGTGGGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTT
...
```

The function accepts the result of the BLAST algorithm.

```python
parse_blast_output(
    input_file="data/example_blast_results.txt",
    output_file="data/blast_output.txt",
)
```

``blast_output.txt`` file contains 37 string with gene names from Description column in input file in alphabetical order (only 10 examples shown):
```python
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

```python
select_genes_from_gbk_to_fasta(
    "data/example_gbk.gbk",
    ["ybcO"],
    n_before=2,
    n_after=2,
    output_fasta="data/TEST2_output_fasta.fasta",
)
```

As a result of the function's operation, ``.fasta`` file is obtained with the following content (only 1 example shown):
```python
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