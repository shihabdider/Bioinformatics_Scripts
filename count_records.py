from Bio.Emboss.Applications import NeedleCommandline
from Bio import SeqIO
import os

def count_records(filename):
    with open('water_outputs\{0}'.format(filename), "r") as handle:
        count = 0
        for line in handle:
            count += 1
        if count != 849:
            print filename

outputs = [name for name in os.listdir("C:\Users\User\Research\GPCR Project\\water_outputs")]

output = 'pairwise_GPCR_water_scores_self_align.txt'
for filename in outputs:
    count_records(filename)
