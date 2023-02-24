# DNA Toolset/Code testing file

from bio_seq import bio_seq
from utilities import readTextFile, writeTextFile, read_FASTA

test_dna = bio_seq()
test_dna.generate_rnd_seq(40, "RNA")


print(test_dna.get_seq_info())
print(f'[Nucleotide Frequency]: {test_dna.count_nuc_frequency()}')
print(f'[DNA/RNA Transcription]: {test_dna.transcription()}\n')

print(f"[DNA String + Complement]:\n\n5' {test_dna.seq} 3'")
print(f"   {''.join(['|' for c in range(len(test_dna.seq))])}")
print(f"3' {test_dna.complement()} 5' [Complement]")

print(f"5' {test_dna.reverse_complement()} 3' [Reverse Complement]\n")
print(f"[GC Content]: {test_dna.gc_content()}%\n")
print(f"[GC Content in Subsection k=20: {test_dna.gc_content_subsec()}\n")
print(f'[Aminoacids Sequence]: {test_dna.translate_seq()}\n')
print(f'[Codon Frequency (L): {test_dna.codon_usage("L")}\n')
print('[Reading Frames]:')
for frame in test_dna.gen_reading_frames():
    print(frame)
print('\n[All proteins in 6 open reading frames]:')
for prot in test_dna.all_proteins_from_open_rfs():
    print(f'{prot}')
