from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

"""
This module provides classes and functions for processing biological sequences, including DNA, RNA, and proteins.
It also includes functionality for filtering FASTQ files based on sequence quality, length, and GC content.

Author: Nika Samusik
"""

DIFF_DNA_RNA = {"T": "U", "t": "u", "U": "T", "u": "t"}
DNA_TO_RNA = {
    "A": "U",
    "T": "A",
    "G": "C",
    "C": "G",
    "a": "u",
    "t": "a",
    "g": "c",
    "c": "g",
}
RNA_TO_RNA = {
    "A": "U",
    "U": "A",
    "G": "C",
    "C": "G",
    "a": "u",
    "u": "a",
    "g": "c",
    "c": "g",
}
DNA_TO_DNA = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
}
RNA_TO_DNA = {
    "A": "T",
    "U": "A",
    "G": "C",
    "C": "G",
    "a": "t",
    "u": "a",
    "g": "c",
    "c": "g",
}
weights = {
            "A": 89,
            "R": 174,
            "N": 132,
            "D": 133,
            "B": 133,
            "C": 121,
            "Q": 146,
            "E": 147,
            "Z": 147,
            "G": 75,
            "H": 155,
            "I": 131,
            "L": 131,
            "K": 146,
            "M": 149,
            "F": 165,
            "P": 115,
            "S": 105,
            "T": 119,
            "W": 204,
            "Y": 181,
            "V": 117,
        }

class BiologicalSequence(ABC):
    """
    An abstract class representing a biological sequence.

    Provides common functionality for handling sequences, such as:
    - Validation of the sequence alphabet.
    - Standardized string representation.
    - Indexing and slicing support.
    """

    def __init__(self, sequence: str, valid_alphabet: set):
        if not sequence:
            raise ValueError("Sequence cannot be empty.")
        self.sequence = sequence.upper()
        self.valid_alphabet = valid_alphabet
        self._validate_sequence()

    def _validate_sequence(self):
        if not set(self.sequence).issubset(self.valid_alphabet):
            raise ValueError(
                f"Invalid symbols: {set(self.sequence) - self.valid_alphabet}"
            )

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.sequence}')"

    @abstractmethod  # переопределяем его дальше
    def get_complement(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    Represents a nucleic acid sequence (DNA or RNA).

    Provides methods for:
    - Complementing the sequence.
    - Reversing the sequence.
    - Obtaining the reverse complement.
    - Calculating the annealing temperature.
    - Checking if the sequence is a valid primer or palindrome.
    """

    complement_map = {}

    def get_complement(self):
        return self.__class__("".join(self.complement_map[n] for n in self.sequence))

    def reverse(self):
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self):
        return self.get_complement().reverse()

    def annealing_temperature(self):
        return 2 * (self.sequence.count("A") + self.sequence.count("T")) + 4 * (
            self.sequence.count("G") + self.sequence.count("C")
        )

    def check_palindrome(self):
        return self.sequence == self.reverse_complement().sequence

    def check_primer(self):
        seq_len = len(self.sequence)
        if seq_len < 16 or seq_len > 30:
            return False
        gc_content = (self.sequence.count("G") + self.sequence.count("C")) / seq_len
        tm = self.annealing_temperature()
        return 0.4 <= gc_content <= 0.6 and tm > 55 and self.sequence[-1] in {"G", "C"}


class DNASequence(NucleicAcidSequence):
    """
    Represents a DNA sequence.

    Provides functionality for transcription into RNA.
    """

    complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def __init__(self, sequence: str):
        super().__init__(sequence, {"A", "T", "G", "C"})

    def transcribe(self):
        return RNASequence("".join(DIFF_DNA_RNA.get(n, n) for n in self.sequence))


class RNASequence(NucleicAcidSequence):
    """
    Represents an RNA sequence.

    Provides functionality for reverse transcription into DNA.
    """

    complement_map = {"A": "U", "U": "A", "G": "C", "C": "G"}

    def __init__(self, sequence: str):
        super().__init__(sequence, {"A", "U", "G", "C"})

    def transcribe(self):
        return DNASequence("".join(DIFF_DNA_RNA.get(n, n) for n in self.sequence))


class AminoAcidSequence(BiologicalSequence):
    """
    Represents an amino acid sequence.

    Provides functionality for computing molecular weight.
    """

    VALID_AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")

    def __init__(self, sequence: str):
        super().__init__(sequence, self.VALID_AA_ALPHABET)

    def molecular_weight(self, weight):
        return sum(weights[aa] for aa in self.sequence)


def filter_fastq(
    input_fastq,
    output_fastq=None,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0,
):
    """
    Filters sequences from a FASTQ file based on GC content, length, and quality.

    Parameters:
    - input_fastq: Path to the input FASTQ file.
    - output_fastq: Optional path to save filtered sequences.
    - gc_bounds: GC content filtering range (default: 0-100).
    - length_bounds: Length filtering range (default: 0-2**32).
    - quality_threshold: Minimum average quality for filtering (default: 0).
    """
    min_gc, max_gc = gc_bounds
    min_len, max_len = length_bounds

    filtered_sequences = []

    with open(input_fastq, "r") as infile:
        records = SeqIO.parse(infile, "fastq")

        for record in records:
            seq_gc = gc_fraction(record.seq) * 100
            seq_len = len(record.seq)
            avg_quality = sum(record.letter_annotations["phred_quality"]) / seq_len

            if (min_gc <= seq_gc <= max_gc and min_len <= seq_len <= max_len and avg_quality >= quality_threshold):
                filtered_sequences.append(record)
    if output_fastq:
        with open(output_fastq, "w") as outfile:
            SeqIO.write(filtered_sequences, outfile, "fastq")
    else:
        return filtered_sequences
