''' Converts SMAP alignment output into a gapped fasta alignment '''

import GPCR_Project_func as GPCR_func

# Process SMAP output
SMAP_dict = GPCR_func.process_SMAP('5cxv_3rze_.txt')

# Load non-gapped fasta files

fasta1 = '3rze.fasta'
aln2 = GPCR_func.records_to_list(fasta1)
fasta2 = '5cxv.fasta'
aln1 = GPCR_func.records_to_list(fasta2)
alns = [aln1[0], aln2[0]]

# Convert SMAP output into fasta alignment 

def SMAP_to_fasta(SMAP_index, alns):
    ''' Inserts gaps into sequences based on index, adjust index as you go '''

    aln1_index = SMAP_index[0]
    aln2_index = SMAP_index[1]
    
    aln1_seq = alns[0].seq.tomutable()
    aln2_seq = alns[1].seq.tomutable()

    diff_aln_1 = 0
    diff_aln_2 = 0
    for i_1, pos_aln_1 in enumerate(aln1_index[:]):
        for i_2, pos_aln_2 in enumerate(aln2_index[:]):
            if i_1 == i_2:
                pos_aln_1 += diff_aln_1
                pos_aln_2 += diff_aln_2
                if pos_aln_1 > pos_aln_2:
                    print 'mismatch 1', pos_aln_1, pos_aln_2
                    diff = pos_aln_1 - pos_aln_2 
                    print 'before', diff_aln_2
                    diff_aln_2 += diff
                    for num in range(diff):
                        aln2_seq.insert(pos_aln_2, '-')
                    print 'after', diff_aln_2

                elif pos_aln_2 > pos_aln_1:
                    print 'mismatch 2', pos_aln_1, pos_aln_2
                    diff = pos_aln_2 - pos_aln_1 
                    print 'before', diff_aln_1
                    diff_aln_1 += diff 
                    for num in range(diff):
                        aln1_seq.insert(pos_aln_1, '-')
                    print 'after', diff_aln_1

                else:
                    print 'no mismatch', pos_aln_1, pos_aln_2

    print aln1_index
    print aln2_index
    print aln1_seq
    print aln2_seq

def gen_SMAP_index(SMAP_dict, aln1, aln2):
    ''' Generates an updated index by accounting for position shifts by
    extracting continous fragments and then using regex pattern matching '''

    import re

    # Extract continous fragments 

    aln_1_re_strings = []
    aln_2_re_strings = []
    
    for i, prot in enumerate(SMAP_dict):
        residue_index = SMAP_dict[prot]
        pos_list = []
        amino_list = []
        for residue in residue_index:
            pos = int(residue[-3:])
            amino = residue[2]

            pos_list.append(pos)
            amino_list.append(amino)

        num_dots_list = [j-k-1 for k, j in zip(pos_list[:-1], pos_list[1:])] 

        regex_list = [None]*(len(amino_list)+len(num_dots_list))
        regex_list[::2] = amino_list
        regex_list[1::2] = num_dots_list
        

        if i == 0:
            re_string = ''
            for elem in regex_list:
                if type(elem) is str:
                    re_string += '({0})'.format(elem)
                elif elem != 0 and elem <= 50:
                    re_string += '.{{{0}}}'.format(elem)
                elif elem > 50:
                    aln_1_re_strings.append([re_string])
                    aln_1_re_strings.append([elem])
                    re_string = ''
            aln_1_re_strings[0].append(re_string)
        elif i == 1:
            re_string = ''
            for elem in regex_list:
                if type(elem) is str:
                    re_string += '({0})'.format(elem)
                elif elem != 0 and elem <= 50:
                    re_string += '.{{{0}}}'.format(elem)
                elif elem > 50:
                    aln_2_re_strings.append([re_string])
                    aln_2_re_strings.append([elem])
                    re_string = ''
            aln_2_re_strings[0].append(re_string)

    # Pattern match the continous fragments and retrieve new positions

    aln1_re_list = make_pattern_from_list(aln_1_re_strings[0])
    aln2_re_list = make_pattern_from_list(aln_2_re_strings[0])

    aln_pos_list = [get_fragment_pos(aln1_re_list, aln1),
            get_fragment_pos(aln2_re_list, aln2)]
        
    return aln_pos_list

def get_fragment_pos(aln_re_list, aln):
    import re

    aln_pos_list = []
    for pattern in aln_re_list:
        m = re.search(pattern, str(aln[0].seq))
        for i in range(1, len(m.groups())):
            aln_pos_list.append([m.group(i), m.start(i)])

    return aln_pos_list

def make_pattern_from_list(pattern_list):
    ''' Makes a list of re objects from a list of pattern strings '''
    import re

    re_list = [re.compile(pattern) for pattern in pattern_list]
    return re_list


def update_SMAP_index(SMAP_index, increment, pos):
    ''' Increments SMAP index '''
    
    new_index = []
    for i, num in enumerate(SMAP_index):
        if i >= pos:
            num += increment
            new_index.append(num)
        else:
            new_index.append(num)

    return new_index
       
new_SMAP_list = gen_SMAP_index(SMAP_dict, aln1, aln2)
print new_SMAP_list
SMAP_index = []
for prot in new_SMAP_list:
    SMAP_index.append([elem[1] for elem in prot])

SMAP_to_fasta(SMAP_index, alns)
