import os


def calculate_gc_bounds(seq) -> float:
    """
    Calculates the GC content of a nucleotide sequence.

    This function computes the percentage of guanine (G) and cytosine (C)
    nucleotides in the given DNA or RNA sequence.

    Parameters
    ----------
    seq : str
        Nucleotide sequence.

    Returns
    -------
    gc_content : float
        The GC content as a percentage of the total sequence.
    """
    count_g = seq.lower().count("g")
    count_c = seq.lower().count("c")
    gc_content = ((count_g + count_c) / len(seq)) * 100
    return gc_content


def calculate_quality_threshold(quality) -> float:
    """
    Calculates the average Q-score of a nucleotide sequence from a FASTQ file.

    This function computes the quality score (Q-score) of a sequence using
    the formula: Unicode value of the character - 33.

    Parameters
    ----------
    quality : str
        A string of characters representing the quality scores of nucleotides
        from a FASTQ file.

    Returns
    -------
    q_score : float
        The average Q-score of the sequence.
    """
    return sum(ord(letter) - 33 for letter in quality) / len(quality)


def make_bounds(inputs) -> tuple:
    """
    Creates bounds from an integer or returns the given bounds as is.

    Parameters
    ----------
    bnds : int or tuple
        The input bounds, either as an integer (upper limit) or a
        tuple (lower limit, upper limit).

    Returns
    -------
    tuple
        A tuple representing the lower and upper bounds.

    Raises
    -------
    ValueError
        If the bounds are not correctly specified.
    """
    if isinstance(inputs, int):
        return (0, inputs)
    elif isinstance(inputs, tuple) and len(inputs) == 2:
        if (
            inputs[0] < 0
            or inputs[1] < 0
            or inputs[0] > inputs[1]
        ):
            raise ValueError("Bounds must be non-negative and\
                            the lower bound must not exceed the upper bound.")
        return inputs
    else:
        raise ValueError("Bounds must be either an integer\
                        or a tuple of two integers.")


def is_bounded(bounds, value) -> bool:
    """
    Checks if a value is within the specified bounds.

    Parameters
    ----------
    bounds : tuple
        A tuple containing the lower and upper bounds.
    value : int or float
        The value to check against the bounds.

    Returns
    -------
    bool
        True if value is within bounds, False otherwise.
    """
    return bounds[0] <= value <= bounds[1]


def read_fastq(input_fastq):
    """
    Reads a FASTQ file and converts it to a dictionary.

    Parameters
    ----------
    input_fastq : str
        The path to the input FASTQ file.

    Returns
    -------
    dict
        A dictionary where the key is the sequence name,
        and the value is a tuple (sequence string, quality string).
    """
    sequences = {}
    with open(input_fastq, 'r') as infile:
        while True:
            header = infile.readline().strip()
            if not header:
                break

            seq = infile.readline().strip()
            infile.readline()
            quality = infile.readline().strip()

            seq_name = header
            sequences[seq_name] = (seq, quality)

    return sequences


def write_fastq(header, seq, quality, output_fastq):
    """
    Writes the filtered sequences to a FASTQ file.

    Parameters
    ----------
    output_fastq : str
        The path to the output FASTQ file.
    filtered_sequences : dict
        A dictionary of filtered sequences where the key is the sequence name
        and the value is a tuple (sequence string, quality string).
    """
    if output_fastq:
        path, filename = os.path.split(output_fastq)
        filtered_path = os.path.join(path, 'filtered')

        os.makedirs(filtered_path, exist_ok=True)
        output_path = os.path.join(filtered_path, filename)

        with open(output_path, 'a') as outfile:
            outfile.write(f"{header}\n{seq}\n{quality}\n")
    else:
        print(f"{header}\n{seq}\n{quality}\n")
