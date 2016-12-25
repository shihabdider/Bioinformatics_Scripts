'''Cleans an alignment file of columns which have a certain percentage of
gaps, then obtains substitution score vectors'''

import numpy
from Bio import pairwise2
from Bio.Seq import MutableSeq
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

records = records_to_list(filename)

test_records = []
test_record_names = ['HRH1', 'ACM1']
for name in test_record_names:
    for record in records:
        if name in record.id:
            test_records.append(record)

gap_list = []

for record in test_records:
    gaps_row = []
    for char in record.seq:
        if char == '-':
            gaps_row.append(0)
        else:
            gaps_row.append(1)
    gap_list.append(gaps_row)

gap_matrix = numpy.asarray(gap_list)

'''Get the columns to remove'''
gaps_limit_ratio = float(1/(len(test_record_names)))

clean_column_index = []
total_column_len = float(len(gap_matrix[:,0]))
for n in range(len(gap_matrix[0,:])):
    total_non_gaps = 0
    for num in gap_matrix[:,n]:
        total_non_gaps += num
    ratio_non_gaps = (total_non_gaps/total_column_len)
    if (ratio_non_gaps)<=gaps_limit_ratio:
        clean_column_index.append(n)

#print clean_column_index, len(records[0].seq)-len(clean_column_index)

'''Remove those columns'''
cleaned_align_index = [pos for pos in range(len(test_records[0].seq)) if pos
        not in clean_column_index]
cleaned_records = []
cleaned_pos_index = []
for record in test_records:
    new_seq = ''
    for pos, char in enumerate(record.seq):
        if pos not in clean_column_index:
            new_seq += char
            cleaned_pos_index.append(pos)
    #break
    #print record.seq, len(record.seq), 'before'
    record.seq = Seq(new_seq)
    #print record.seq, len(record.seq), 'after'
    cleaned_records.append(record)

#print cleaned_pos_index
print cleaned_records
#print 'FInished cleaning records, calculating vectors...'

def write_pairwise_scores(data, filename):
    import csv
    with open(filename, 'ab') as output:
        csv_obj = csv.writer(output, delimiter=',')
        csv_obj.writerows([data])

with open('HRH1_ACM1_gross.fasta', 'a') as handle:
    for record in cleaned_records:
        SeqIO.write(record, handle, 'fasta')

#write_pairwise_scores(cleaned_pos_index, 'cleaned_pos_index.csv')


#matrix = matlist.blosum62
#matrix['-', '-'] = 0
#for char in 'ARNDCQEGHILKMFPSTWYV':
#    matrix[char, '-'] = -10
#    matrix['-', char] = -10
#
#from itertools import combinations
#record_pairs = list(combinations([record for record in cleaned_records], 2))
##record_pairs = []
##for record1 in records:
##    for record2 in records:
##        if record1.id == record2.id:
##            record_pairs.append([record1, record2])
#
#for record_pair in record_pairs:
#    score_vector = []
#    for n, sequence in enumerate(record_pair[0].seq):
#        align = pairwise2.align.localdx(record_pair[0].seq[n],
#                record_pair[1].seq[n], matrix, score_only=True)
#        print len(score_vector)
#        score_vector.append(align)
#        #score_vector.append([cleaned_align_index[n], align])
#    write_pairwise_scores(score_vector,
#            'pairwise_aa_substitution_vectors_new.csv' )
#    #write_pairwise_scores([record1.id, record2.id, score_vector],
#    #'pairwise_aa_substitution_vectors.csv' )
#
