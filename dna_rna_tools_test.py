from bio_tools import run_dna_rna_tools as rdrt


def test_transcribe():
    assert rdrt("ATG", "transcribe") == "AUG"
    assert rdrt("AGt", "transcribe") == "AGu"


def test_reverse():
    assert rdrt("ATG", "reverse") == "GTA"
    assert rdrt("cUG", "reverse") == "GUc"


def test_complement():
    assert rdrt("AtG", "complement") == "TaC"
    assert rdrt("CUG", "complement") == "GAC"


def test_reverse_complement():
    assert rdrt("ATg", "reverse_complement") == "cAT"
    assert rdrt("CUG", "reverse_complement") == "CAG"


def test_multiple_args():
    assert rdrt("ATG", "aT", "reverse") == ["GTA", "Ta"]
    assert rdrt("ttG", "AT", "ATc", "complement") == ["aaC", "TA", "TAg"]


def test_palindrom():
    assert rdrt("ACCGCGGT", "ACCGCGGT", "is_palindrome") == [True, True]
    assert rdrt("ACCGCGGT", "AT", "ATc", "is_palindrome") == [True, True, False]


def test_primer():
    assert rdrt("GTTGTAAAACGACGGCCAGTGGGG", "AGCGGATAACAATTTCACACAGGAGGGGC", "is_primer") == [True, True]
    assert rdrt("GAaTTgAATTc", "ACCGCGGT", "AT", "ATc", "is_primer") == [False, False, False, False]


def test_annealing_temperature():
    assert rdrt("ATGC", "atGc", "annealing_temperature") == [12, 12]
    assert rdrt("atGc", "annealing_temperature") == 12