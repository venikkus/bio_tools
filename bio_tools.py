from modules import dna_rna_tools_modules as drtm
from modules import filter_fastq_modules as ffm

"""
The module provides functions for processing
DNA/RNA sequences and filtering FASTQ data.

Functions:
- run_dna_rna_tools: Performs operations on DNA/RNA sequences.
- filter_fastq: Filters FASTQ files by specified
parameters of quality, length and GC content.

Author: Nika Samusik
"""


def run_dna_rna_tools(*args) -> list | str:
    """
    Performs an operation on DNA or RNA sequences
from the drtm module

    Parameters
    ----------
    *args : tuple, str
        DNA/RNA sequences and an element with the operation name.
        The operation must be passed as the last argument.

    Returns
    -------
    answer : list, str
        Results or result of executing functions from the module
        drtm. Invalid input warning
        if the submitted data does not contain an operation or
        submits an invalid operation

    Raises
    -------
    ValueError
        if the sequences are not DNA/RNA or
    KeyError
        if the operation is not supported
    """
    if len(args) == 1:
        raise ValueError("There is no operation or sequence")

    operation = args[-1]
    args = args[:-1]

    # check if it is a DNA or RNA at all
    if not drtm.is_nucleotide(*args):
        raise ValueError("This is not a DNA/RNA sequence at all")

    try:
        answer = [getattr(drtm,
                          operation)(arg) for arg in args]
    except KeyError:
        raise KeyError(f"Operation {operation} is not supported.")
    else:
        return answer[0] if len(answer) == 1 else answer


def filter_fastq(
    input_fastq,
    output_fastq=None,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0
):
    """
    Filters sequences according to quality, length and GC content parameters.

    Parameters
    ----------
    fastq_file : dict
        input data: dictionary, where key is the sequence name,
        value is a tuple (sequence string, quality string).
    gc_bounds : int, tuple, default: 0, 100
        GC content thresholds in the sequence.
        You can set an upper threshold and not set a lower one
    length_bounds : int, tuple, default: 0, 2**32
        sequence length thresholds.
        You can set an upper threshold and not set a lower one
    quality_threshold : int, default: 0
        Minimum quality threshold for filtering sequences

    Returns
    -------
    dict
        Filtered FASTQ data

    Raises
    -------
    ValueError
        if the GC composition or length values ​​are outside the limits.
    """
    gc_bounds = ffm.make_bounds(gc_bounds)
    length_bounds = ffm.make_bounds(length_bounds)

    if output_fastq:
        outfile = open(output_fastq, 'w')
    else:
        outfile = None

    with open(input_fastq, 'r') as infile:
        while True:
            header = infile.readline().strip()
            if not header:
                break

            seq = infile.readline().strip()
            infile.readline()  # skip the quality header.
            quality = infile.readline().strip()

            gc_content = ffm.calculate_gc_bounds(seq)
            seq_len = len(seq)
            q_score = ffm.calculate_quality_threshold(quality)

            if (
                q_score >= quality_threshold
                and ffm.is_bounded(gc_bounds, gc_content)
                and ffm.is_bounded(length_bounds, seq_len)
            ):
                ffm.write_fastq(header, seq, quality, output_fastq)
    if outfile:
        outfile.close()


filter_fastq(
    input_fastq='data/example_fastq.fastq',
    output_fastq='data/output_fastq.txt',
    gc_bounds=(40, 100),
    length_bounds=(15, 2**32),
    quality_threshold=33
)
