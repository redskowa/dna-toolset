from collections import Counter
import random

from bio_structs import NUCLEOTIDE_BASE, DNA_Codons, RNA_Codons

class bio_seq:
    """DNA sequence class. Default value: ATCG, DNA, No Label"""

    def __init__(self, seq="ATCG", seq_type="DNA", label="No Label"):
        """Sequence initialization, validation"""
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data does not seem to be correct {self.seq_type} sequence"

    # DNA Toolkit functions
    def __validate(self):
        """Check the sequence to make sure it is a valid DNA string"""
        return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)

    def get_seq_biotype(self):
        """Returns sequence type"""
        return self.seq_type

    def get_seq_info(self):
        """Returns 4 strings. Full sequence information"""
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}"

    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        """Generate a random DNA or RNA sequence, provided the length"""
        seq = ''.join([random.choice(NUCLEOTIDE_BASE[seq_type])
                       for x in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")

    # === Nucleotide Frequency Counter === #
    def count_nuc_frequency(self):
        """Counting nucleotide frequency of a given sequence."""
        return dict(Counter(self.seq))

    # === Transcription === #
    def transcription(self):
        """DNA -> RNA Transcription, Replacing Thymine with Uracil"""
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        return "Not a DNA sequence"

    # === Complement === #
    def complement(self):
        """Swapping adenine with thymine and guanine with cytosine. 3' to 5' order"""
        #return ''.join([DNA_Complement[nuc] for nuc in self.seq])
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG', 'TAGC')
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)

    # === Reverse Complement === #
    def reverse_complement(self):
        '''Swapping adenine with thymine and guanine with cytosine. Then reversing generated string. 5' to 3' order'''
        # return ''.join([DNA_Complement[nuc] for nuc in self.seq])[::-1]

        #pythonic solution using built in .maketrans function
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG', 'TAGC')
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)[::-1]

    ## ==== GC Content ==== ##
    def gc_content(self):
        '''GC content ration in a DNA/RNA sequence'''
        return round((self.seq.count('C') + self.seq.count('G')) / len(self.seq) * 100, 3)

    def gc_content_subsec(self, k=20):
        '''GC Content in a DNA/RNA sub-sequence length k. k = 20 by default'''
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i+k]
            res.append(
                round((subseq.count('C') + subseq.count('G')) / len(subseq) * 100, 3)
                )
        return res

    def translate_seq(self, init_pos=0):
        '''Translates a DNA sequence into an amino acid sequence'''
        if self.seq_type == "DNA":
            return [DNA_Codons[self.seq[pos:pos+3]] for pos in range(init_pos, len(self.seq) -2, 3)]
        elif self.seq_type == "RNA":
            return [RNA_Codons[self.seq[pos:pos+3]] for pos in range(init_pos, len(self.seq) -2, 3)]

    def codon_usage(self, aminoacid):
        '''Provides the frequency of each codon encoding a given aminoacid in a DNA sequence'''
        tmpList = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) -2, 3):
                if DNA_Codons[self.seq[i:i+3]] == aminoacid:
                    tmpList.append(self.seq[i:i+3])
        elif self.seq_type == "RNA":
            for i in range(0, len(self.seq) -2, 3):
                if RNA_Codons[self.seq[i:i+3]] == aminoacid:
                    tmpList.append(self.seq[i:i+3])
        
        freqDict = dict(Counter(tmpList))
        totalWight = sum(freqDict.values())
        for self.seq in freqDict:
            freqDict[self.seq] = round(freqDict[self.seq] / totalWight, 2)
        return freqDict

    def gen_reading_frames(self):
        """Generate the six reading frames of a DNA sequence, including the reverse complement"""
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tmp_seq = bio_seq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        del tmp_seq
        return frames

    def proteins_from_rf(self, aa_seq):
        '''Compute all possible proteins in an aminoacid seq and return a list of possible proteins'''
        current_prot = []
        proteins = []
        for aa in aa_seq:
            if aa == "_":
                # STOP accumulating amino acids if _ - STOP was found
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                # START accumulating amino acids if M - START was found
                if aa == "M":
                    current_prot.append("")
                for i in range (len(current_prot)):
                    current_prot[i] += aa
        return proteins

    def all_proteins_from_open_rfs(self, startReadPos=0, endReadPos=0, ordered=False):
        '''Compute all possible proteins for all open reading frames'''
        '''Protein search DB: https://www.ncbi.nlm.nih.gov/nucore/NM_001185097.2'''
        '''API can be used to pull protein info'''
        if endReadPos > startReadPos:
            tmp_seq = bio_seq(self.seq[startReadPos: endReadPos], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
            del tmp_seq
        else:
            rfs = self.gen_reading_frames()

        result = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                result.append(p)
        
        if ordered:
            return sorted(result, key=len, reverse=True)
        
        return result