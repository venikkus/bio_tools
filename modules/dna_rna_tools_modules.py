DIFF_DNA_RNA = {
    "T": "U",
    "t": "u",
    "U": "T",
    "u": "t"
}
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


# transcribe — вернуть транскрибированную последовательность
def transcribe(arg) -> str:
    """
    Transcribes the given DNA or RNA sequence.

    This function replaces nucleotides in the sequence based on transcription
    rules. For example, DNA adenine (A) is transcribed to RNA uracil (U), and
    vice versa.

    Parameters
    ----------
    arg : str
        A nucleotide sequence: either DNA or RNA.

    Returns
    -------
    result : str
        The transcribed RNA or DNA nucleotide sequence.
    """
    result = ""

    for letter in arg:
        if letter in DIFF_DNA_RNA:
            result += DIFF_DNA_RNA[letter]
        else:
            result += letter
    return result


# reverse — вернуть развёрнутую последовательность
def reverse(arg) -> str:
    """
    Reverses the given DNA or RNA sequence.

    This function returns the reversed version of the input nucleotide
    sequence.

    Parameters
    ----------
    arg : str
        A nucleotide sequence: either DNA or RNA.

    Returns
    -------
    result : str
        The reversed nucleotide sequence.
    """
    return arg[::-1]


# complement — вернуть комплементарную последовательность
def complement(arg) -> str:
    """
    Returns the complementary sequence for the given DNA or RNA sequence.

    This function computes the complementary nucleotide sequence by
    replacing each nucleotide with its pair: A with T/U,
    C with G, and vice versa.

    Parameters
    ----------
    arg : str
        A nucleotide sequence: either DNA or RNA.

    Returns
    -------
    result : str
        The complementary nucleotide sequence.
    """
    result = ""
    for letter in arg:
        if letter in DNA_TO_DNA:
            result += DNA_TO_DNA[letter]
        else:
            result += RNA_TO_RNA[letter]
    return result


# reverse_complement — вернуть обратную комплементарную последовательность
def reverse_complement(arg) -> str:
    """
    Returns the reverse complement of the given DNA or RNA sequence.

    This function computes both the complement and the reverse
    of the input sequence.

    Parameters
    ----------
    arg : str
        A nucleotide sequence: either DNA or RNA.

    Returns
    -------
    result : str
        The reverse complement of the nucleotide sequence.
    """
    arg = complement(arg)
    return reverse(arg)


# is_palindrome — является ли поданная последовательность биопалиндромом
def is_palindrome(arg) -> bool:
    """
    Checks whether the given sequence is a biological palindrome.

    A biological palindrome is a sequence that is equal to
    it's reverse complement.

    Parameters
    ----------
    arg : str
        A nucleotide sequence.

    Returns
    -------
    bool
        True if the sequence is a biological palindrome, False otherwise.
    """
    reverse = ""
    for letter in arg.lower():
        if letter in DNA_TO_DNA:
            reverse += DNA_TO_DNA[letter]
        else:
            reverse += RNA_TO_RNA[letter]

    return arg.lower() == reverse[::-1].lower()


# annealing_temperature — подсчёт температуры отжига праймера
def annealing_temperature(arg) -> float:
    """
    Calculates the melting temperature (Tm) of the given sequence.

    The Tm is calculated using the formula:
    Tm (°C) = 2 * (A + T) + 4 * (G + C)

    Parameters
    ----------
    arg : str
        A nucleotide sequence.

    Returns
    -------
    result : float
        The melting temperature of the sequence.
    """
    arg = arg.lower()
    count_a = arg.count("a")
    count_t = arg.count("t")
    count_g = arg.count("g")
    count_c = arg.count("c")

    # Tm (°C ) = 2 х (A+T) + 4 х (G+C)
    result = (2 * (count_a + count_t)) + (4 * (count_g + count_c))
    return result


# is_primer — является ли поданная последовательность праймером
def is_primer(arg) -> bool:
    """
    Checks if the given sequence meets primer criteria.

    A primer must have a GC content between 40% and 60%
    and a melting temperature (Tm) > 55°C.

    Parameters
    ----------
    arg : str
        A nucleotide sequence.

    Returns
    -------
    bool
        True if the sequence meets primer criteria, False otherwise.
    """
    arg = arg.lower()

    if len(arg) < 16 or len(arg) > 30 or arg[-1].lower() not in ["g", "c"]:
        return False
    else:
        count_g = arg.count("g")
        count_c = arg.count("c")

        temp = annealing_temperature(arg)
        gc = (count_g + count_c) / len(arg)

    return 0.4 < gc < 0.6 or 55 < temp


# is_nucleotide — является ли поданная последовательность нуклеотидной
def is_nucleotide(*args) -> bool:
    """
    Checks if the given sequence(s) represent valid DNA or RNA.

    A valid nucleotide sequence contains only A, T, G, C (for DNA)
    or A, U, G, C (for RNA), and it cannot contain both U and T
    in the same sequence.

    Parameters
    ----------
    *args : str
        A list of nucleotide sequences.

    Returns
    -------
    bool
        True if all provided sequences are valid DNA or RNA, False otherwise.
    """
    for arg in args:
        if "u" in arg.lower() and "t" in arg.lower():
            return False
        for letter in arg:
            if letter.lower() not in "atgcu":
                return False
    return True
