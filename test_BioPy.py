import sys

import numpy as np

from Bio import SeqIO
from Bio import Seq
from array import array
from pprint import pprint


# fasta_rec_list = list(SeqIO.parse("data\example.fasta", "fasta"))
fasta_list = list(SeqIO.parse("data\example.fasta", "fasta"))
# pprint(fasta_list)
# print(fasta_list[0].id)
full_str = ''.join(seq_rec.seq._data for seq_rec in fasta_list)
print(full_str)
