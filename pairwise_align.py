from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
import os

#seq_fnames = [name for name in os.listdir("C:\Users\User\Research\GPCR Project\GPCRs Binding Fastas")]
#
#seq_fname2 = 'C:\Users\User\Research\GPCR Project\\binding_site_aln.fasta'
## needle_fname = 'C:\Users\User\Research\GPCR Project\pairwise_binding_output.emboss'
#
#count = 0
#for seq_fname in seq_fnames:
#    seq_fname1 = 'C:\Users\User\Research\GPCR Project\GPCRs Binding Fastas\{0}'.format(seq_fname)
#    needle_fname = 'C:\Users\User\Research\GPCR Project\\bs_outputs\{0}'.format(seq_fname)
#    needle_cli = NeedleCommandline(asequence=seq_fname1, \
#                               bsequence=seq_fname2, \
#                               gapopen=10, \
#                               gapextend=0.5, \
#                               aformat='score',
#                               outfile=needle_fname)
#
#    needle_cli() 
#    count += 1
#    print 'Parsed {0}'.format(seq_fname), count
#"""This generates the needle file"""
#"""That parses the needle file, aln[0] and aln[1] contain the aligned
#first and second sequence in the usual format (e.g. - for a gap)"""
#

def parse_emboss_scores(filename, output):
    with open(output, 'a') as output:
        with open('water_outputs\{0}'.format(filename), "r") as handle:
            for line in handle:
                words = line.split()
                if len(words)>2:
                    if words[0] != words[1]:
                        output.write(words[0]+'\t')
                        output.write(words[1]+'\t')
                        output.write(words[3]+'\t')
                        output.write('\n')

outputs = [name for name in os.listdir("C:\Users\User\Research\GPCR Project\\water_outputs")]

output = 'pairwise_GPCR_water_scores.txt'
for filename in outputs:
    parse_emboss_scores(filename, output)
