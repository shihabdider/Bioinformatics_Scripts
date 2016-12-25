'''Extracts positions using gap/non-gap frames of reference'''

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import GPCR_Project_func as GPCR_func


def extract_tm_region(filename):
    with open(filename, 'r') as tm_region:
        residue_pos = []
        for i, line in enumerate(tm_region):
            if i >= 2:
                line_list = line.split()
                for num in range(int(line_list[0]), int(line_list[1])+1):
                    residue_pos.append(num)

        return list(set(residue_pos))


def extract_pos(gaps, records, residue_index_dict, start_pos=1):
    '''Extracts the specified residue index positions from a different frame of
    reference (if gap, then from non-gap and vice versa)

    Output
    - returns a dictioanry whose keys = record.ids and values = [[desired pos,
      given pos, residue],]
    
    Arguments
    - gaps: Whether the sequence to be extracted from has gaps or not
    - records: The sequence records which will be used
    - residue_index_dict: A dictionary containing the id of the protein (key)
      and the residues to be extracted (value)'''

    pos_dict = {}
    for entry in residue_index_dict:
        if gaps:
            pos_list = extract_pos_gaps(records, residue_index_dict[entry],
                    entry, start_pos)
            pos_dict[entry] = pos_list

        else:
            pos_list = extract_pos_no_gaps(records, residue_index_dict[entry],
                    entry, start_pos)
            pos_dict[entry] = pos_list

    return pos_dict

def extract_pos_gaps(records, residue_index, template_id, start_pos):
    '''Extracts residues given from a gapped sequence, given non-gapped residue
    index
    
    Output
    - Returns a list of positions corresponding to the residues extracted in
      the form: [record.id, [[residue, gapped position, non-gapped position],]]

    See extract_pos for arguments'''
    
    pos_index = []
    for record in records:
        if template_id in record.id:
            print record.id, 'Extracting specified residues...'
            gapped_count = -1
            seq_count = start_pos
            
            site_list = []
            for char in record.seq:
                gapped_count += 1
                if char != '-':
                    if seq_count in residue_index:
                        #print ('''Residue: {0}, Gapped Position: {1}, Ungapped
                        #Position: {2}'''.format(record.seq[gapped_count],
                        #    gapped_count, seq_count))
                        site_list.append([gapped_count, seq_count, char])
                    seq_count += 1

            pos_index = [record.id, site_list]
            return pos_index

def extract_pos_no_gaps(records, residue_index, template_id, start_pos):
    '''Extract non-gapped positions given gapped positions
    
    See extract_pos_gaps for arguments and output (just flip gapped with
    non-gapped)'''

    pos_index = []
    for record in records:
        if template_id in record.id:
            print record.id, 'Extracting specified residues...'
            gapped_count = -1
            seq_count = start_pos
            
            site_list = []
            for char in record.seq:
                gapped_count += 1
                if char != '-':
                    if gapped_count in residue_index:
                        #print '''Residue: {0}, Ungapped Position: {1}, Gapped
                        #Position: {2}'''.format(record.seq[gapped_count],
                        #        seq_count, gapped_count)
                        site_list.append([seq_count, gapped_count, char])
                    seq_count += 1

            pos_index = [record.id, site_list]
            return pos_index

def get_pos_seq_dict(seq_index_dict, prot_id):
    ''' Get the positions from the position index return by the extract_pos
    function for a specified protien'''

    for key in seq_index_dict:
        if prot_id in key:
            residue_index = seq_index_dict[key][1]
            return [pos[0] for pos in residue_index]    

# Load records in fasta to a list
filename = 'HRH1_ACM1_gross.fasta'
filename_fatcat = 'HRH1_ACM1_FATCAT.fasta'
records = GPCR_func.records_to_list(filename)
records_fatcat = GPCR_func.records_to_list(filename_fatcat)

training_prots_records = GPCR_func.records_to_list(
    'top_ranked_residues_training_set.fasta')

# Specify the positions to be extracted
tm_region_HRH1 = extract_tm_region('tm_region_HRH1.txt')
tm_region_ACM1 = extract_tm_region('tm_region_ACM1.txt')


residue_index_dict = {
        'HRH1': tm_region_HRH1,
        'ACM1': tm_region_ACM1,}

# Extract those positions depending on frame of reference (gaps or no gaps)

# HRH1 #
########

print 'HRH1 (3rze) TM site\n'

gross_seq_index_dict = extract_pos(True, records, residue_index_dict,
        start_pos=1)

fatcat_seq_index_dict = extract_pos(True, records_fatcat, residue_index_dict,
        start_pos=28)

print '\nHRH1 gross seq dict: ', gross_seq_index_dict['HRH1'], '\n'
print 'HRH1 fatcat seq dict: ', fatcat_seq_index_dict['HRH1'], '\n'

key_id = []
gross_gapped_HRH1 = get_pos_seq_dict(gross_seq_index_dict, 'HRH1')
fatcat_gapped_HRH1 = get_pos_seq_dict(fatcat_seq_index_dict, 'HRH1')

#print gross_gapped_HRH1
#print fatcat_gapped_HRH1

pairs_gross_HRH1 = []
for pos in gross_gapped_HRH1:
    pairs_gross_HRH1.append([records[0].seq[pos],
        records[1].seq[pos]])

pairs_fatcat_HRH1 = []
for pos in fatcat_gapped_HRH1:
    pairs_fatcat_HRH1.append([records_fatcat[0].seq[pos],
        records_fatcat[1].seq[pos]])

pairs_gross_HRH1 = pairs_gross_HRH1[12:]

#print 'HRH1 binding site pairs for gross align: \n', pairs_gross_HRH1, '\n'
#print 'HRH1 binding site pairs for fatcat align: \n', pairs_fatcat_HRH1, '\n'

#pairs_gross_HRH1 = pairs_gross_HRH1[:-1]
print len(pairs_gross_HRH1)

pairs_gross_HRH1_final = []
for i, pair in enumerate(pairs_gross_HRH1):
    if pairs_gross_HRH1[i][0] == pairs_fatcat_HRH1[i][0]:
        pairs_gross_HRH1_final.append(pair)
    #else:
    #    print 'Mismatch found! ', i, pair, pairs_fatcat_HRH1[i],
    #    print gross_seq_index_dict['HRH1'][1][i+12], '\n'
    #    #print fatcat_seq_index_dict['HRH1'][1][i]

pairs_fatcat_HRH1_final = []
for i, pair in enumerate(pairs_fatcat_HRH1):
    if pairs_gross_HRH1[i][0] == pairs_fatcat_HRH1[i][0]:
        pairs_fatcat_HRH1_final.append(pair)

print len(pairs_gross_HRH1_final), 
num_matches = 0.0
for i, pair in enumerate(pairs_gross_HRH1_final):
    if pairs_gross_HRH1_final[i] == pairs_fatcat_HRH1_final[i]:
        #print 'Match found!'
        #print pairs_gross_HRH1[i], pairs_fatcat_HRH1[i]
        num_matches += 1

print 'HRH1 tm site pairs for gross align: \n', pairs_gross_HRH1_final, '\n'
print 'HRH1 tm site pairs for fatcat align: \n', pairs_fatcat_HRH1_final, '\n'

print 'HRH1 pair similarity ratio: ', num_matches/len(pairs_gross_HRH1_final), '\n'
#
# ACM1 #
########
print 'ACM1 (5cvx) TM site \n'

gross_seq_index_dict = extract_pos(True, records, residue_index_dict,
        start_pos=1)

fatcat_seq_index_dict = extract_pos(True, records_fatcat, residue_index_dict,
        start_pos=26)

print '\nACM1 gross seq dict: ', gross_seq_index_dict['ACM1'], '\n'
print 'ACM1 fatcat seq dict: ', fatcat_seq_index_dict['ACM1'], '\n'

gross_gapped_ACM1 = get_pos_seq_dict(gross_seq_index_dict, 'ACM1')
fatcat_gapped_ACM1 = get_pos_seq_dict(fatcat_seq_index_dict, 'ACM1')

#print gross_gapped_ACM1
#print fatcat_gapped_ACM1

pairs_gross_ACM1 = []
for pos in gross_gapped_ACM1:
    pairs_gross_ACM1.append([records[1].seq[pos],
        records[0].seq[pos]])

pairs_fatcat_ACM1 = []
for pos in fatcat_gapped_ACM1:
    pairs_fatcat_ACM1.append([records_fatcat[1].seq[pos],
        records_fatcat[0].seq[pos]])


pairs_gross_ACM1 = pairs_gross_ACM1[10:]
print len(pairs_gross_ACM1), len(pairs_fatcat_ACM1)

pairs_gross_ACM1_final = []
for i, pair in enumerate(pairs_gross_ACM1):
    if pairs_gross_ACM1[i][0] == pairs_fatcat_ACM1[i][0]:
        pairs_gross_ACM1_final.append(pair)
    #else:
    #    print 'Mismatch found! ', i, pair, pairs_fatcat_ACM1[i],
    #    print gross_seq_index_dict['ACM1'][1][i+10], '\n'
    #    #print fatcat_seq_index_dict['HRH1'][1][i]

pairs_fatcat_ACM1_final = []
for i, pair in enumerate(pairs_fatcat_ACM1):
    if pairs_gross_ACM1[i][0] == pairs_fatcat_ACM1[i][0]:
        pairs_fatcat_ACM1_final.append(pair)

#print tm_region_ACM1

print len(pairs_gross_ACM1_final), len(pairs_fatcat_ACM1_final)
num_matches = 0.0
for i, pair in enumerate(pairs_gross_ACM1_final):
    if pairs_gross_ACM1_final[i] == pairs_fatcat_ACM1_final[i]:
        #print 'Match found!'
        #print pairs_gross_ACM1[i], pairs_fatcat_ACM1[i]
        num_matches += 1

print 'ACM1 binding site pairs for gross align: \n', pairs_gross_ACM1, '\n'
print 'ACM1 binding site pairs for fatcat align: \n', pairs_fatcat_ACM1, '\n'

print 'ACM1 pair similarity ratio: ', num_matches/len(pairs_gross_ACM1_final), '\n'

#gapped_index_dict = {
#        'ACM1': gapped_index}
#seq_gapped_index_dict = extract_pos(False, records, gapped_index_dict,
#        start_pos=1)
#
#print seq_gapped_index_dict
#
#
#
