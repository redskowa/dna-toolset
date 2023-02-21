# DNA Toolkit file
import collections
Nucleotides = ["A", "C", "G", "T"]

## Validates the Sequence ##
# Check the sequence to make sure it is a DNA string
def validateSeq (dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq

## Nucleotide Frequency Counter ##
def countNucFrequency(seq):
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict

# Alternate way to code the counter using collections
# def countNucFrequency(seq):
#     return dict(collections.Counter(seq))
