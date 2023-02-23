from bio_structs import *

class bio_seq:
    '''DNA sequence class. Default value: ATCG, DNA, No Label'''

    def __init__(self, seq="ATCG", seq_type="DNA", label="No Label"):
        '''Sequence initialization, validation'''
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data does not seem to be correct {self.seq_type} sequence"
    
    def validate(self):
        '''Check the sequence to make sure it is a valid DNA string'''
        return set(DNA_Nucleotides).issuperset(self.seq)


    def get_seq_biotype(self):
        '''Returns sequence type'''
        return self.seq_type

    
    def get_seq_info(self):
        '''Returns 4 strings. Full sequence information'''
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}"
