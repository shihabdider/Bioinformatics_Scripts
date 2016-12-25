'''Gets the pairwise substitution scores for ALL amino acids as a vector'''

from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO
import os
from Bio.SubsMat import MatrixInfo as matlist

filename = 'gross-alignment.fasta'
def records_to_list(filename):
    records = []
    with open(filename, 'rU') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            records.append(record)
    return records

def write_pairwise_scores(data, filename):
    with open(filename, 'a') as output:
        for elem in data:
            output.write(str(elem))
            output.write('\t')
        output.write('\n')

records = records_to_list(filename)

matrix = matlist.blosum62
matrix['-', '-'] = 0
for char in 'ARNDCQEGHILKMFPSTWYV':
    matrix[char, '-'] = -10
    matrix['-', char] = -10

for record1 in records:
    for record2 in records:
        if record1.id != record2.id:
            score_vector = []
            for n in range(len(record1.seq)-1):
                align = pairwise2.align.localdx(record1.seq[n], record2.seq[n],
                        matrix, score_only=True)
                if align != 0:
                    score_vector.append(align)
            write_pairwise_scores([record1.id, record2.id, score_vector],
            'new_pairwise_scores.txt' )


