# DNA Toolset/Code testing file

from bio_seq import bio_seq

test_dna = bio_seq("ATTGGCAGGTACGTAGGC", "DNA", "Test Label")

print(test_dna.get_seq_info())
print(test_dna.get_seq_biotype())

