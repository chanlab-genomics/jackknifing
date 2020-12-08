from Bio import SeqIO
from Bio import Seq
from array import array
from pprint import pprint


fasta_dict = SeqIO.to_dict(SeqIO.parse("fasta_files\example.fasta", "fasta"))
pprint(fasta_dict)

# first_seq = fasta_dict['Snec_CCMP2469.scaffold1'].seq
first_seq_array = array('u', fasta_dict['Snec_CCMP2469.scaffold1'].seq)
del first_seq_array[-100:]
print(len(first_seq_array))
first_seq_array = Seq.Seq(str(first_seq_array))
