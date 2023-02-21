# DNA Toolset/Code testing file
from DNAToolkit import *
import random

randDNAstr = ''.join([random.choice(Nucleotides)
                    for nuc in range(50)])

DNAStr = validateSeq(randDNAstr)

print(DNAStr)
print(countNucFrequency(DNAStr))
