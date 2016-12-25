'''Gets the similarity score based on strict alignment (meaning no sliding and
then taking the minimum) of the binding scores'''

from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO
import os
from Bio.SubsMat import MatrixInfo as matlist

filename = 'binding_site_aln.fasta'
def records_to_list(filename):
    records = []
    with open(filename, 'rU') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            records.append(record)
    return records

def write_pairwise_scores(data, filename):
    with open(filename, 'a') as output:
        output.write(str(data[0]))
        output.write('\t')
        output.write(str(data[1]))
        output.write('\t')
        output.write(str(data[2]))
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
        if record1.id == record2.id:
            total_score = 0
            for n in range(len(record1.seq)-1):
                align = pairwise2.align.localdx(record1.seq[n], record2.seq[n],
                        matrix, score_only=True)
                total_score += align
            write_pairwise_scores([record1.id, record2.id, total_score],
            'pairwise_bs_strict_self_align_output.txt' )


