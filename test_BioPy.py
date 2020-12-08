import sys

import numpy as np

from Bio import SeqIO
from Bio import Seq
from array import array
from pprint import pprint


# fasta_rec_list = list(SeqIO.parse("data\example.fasta", "fasta"))
fasta_dict = SeqIO.to_dict(SeqIO.parse("data\example.fasta", "fasta"))

# first_seq_array = array('u', fasta_dict['Snec_CCMP2469.scaffold2'].seq)
first_seq_array = fasta_dict['Snec_CCMP2469.scaffold2'].seq._data
first_seq_array = first_seq_array[:100] + first_seq_array[300:]
# print((fasta_dict['Snec_CCMP2469.scaffold2'].seq.__dict__))
# first_seq_array = Seq.Seq(str(first_seq_array))
fasta_dict['Snec_CCMP2469.scaffold2'].seq._data = first_seq_array

# Write the new output file
# SeqIO.write(fasta_rec_list, sys.stdout, 'fasta')
# print(fasta_dict)

record_list = []

for rec_id, rec_val in zip(fasta_dict.keys(), fasta_dict.values()):
    rec_val.id = rec_id
    record_list.append(rec_val)

# SeqIO.write(record_list, sys.stdout, 'fasta')
