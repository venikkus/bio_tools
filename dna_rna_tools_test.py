from bio_tools import DNASequence, RNASequence
import pytest


def test_transcribe():
    assert str(DNASequence("ATG").transcribe()) == "AUG"
    assert str(DNASequence("AGt").transcribe()) == "AGU"


def test_reverse():
    assert str(DNASequence("ATG").reverse()) == "GTA"
    assert str(RNASequence("cUG").reverse()) == "GUC"


def test_complement():
    assert str(DNASequence("AtG").get_complement()) == "TAC"
    assert str(RNASequence("CUG").get_complement()) == "GAC"


def test_reverse_complement():
    assert str(DNASequence("ATg").reverse_complement()) == "CAT"
    assert str(RNASequence("CUG").reverse_complement()) == "CAG"


def test_multiple_args():
    assert [str(DNASequence(seq).reverse()) for seq in ["ATG", "aT"]] == ["GTA", "TA"]
    assert [str(DNASequence(seq).get_complement()) for seq in ["ttG", "AT", "ATc"]] == ["AAC", "TA", "TAG"]


def test_palindrom():
    assert [str(DNASequence(seq).check_palindrome()) for seq in ["ACCGCGGT", "ACCGCGGT"]] == ['True', 'True']
    assert [str(DNASequence(seq).check_palindrome()) for seq in ["ACCGCGGT", "AT", "ATc"]] == ['True', 'True', 'False']


def test_primer():
    assert [str(DNASequence(seq).check_primer()) for seq in ["GTTGTAAAACGACGGCCAGTGGGG", 
                                                             "AGCGGATAACAATTTCACACAGGAGGGGC"]] == ['True', 'True']
    assert [str(DNASequence(seq).check_primer()) for seq in ["GAaTTgAATTc", 
                                                             "ACCGCGGT", "AT", "ATc"]] == ['False', 'False', 'False', 'False']


def test_annealing_temperature():
    assert DNASequence("ATGC").annealing_temperature() == 12
    assert DNASequence("atGc").annealing_temperature() == 12


def test_invalid_dna():
    DNASequence("ATGX")


def test_empty_string():
    DNASequence("")
