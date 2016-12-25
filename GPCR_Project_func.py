from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def process_SMAP(SMAP_output):
    ''' Returns a dictionary with name of the proteins and their alignments '''
    import re
    
    re_dict = {
        'templ_prot_re' : re.compile('Template=.+_'), 
        'query_prot_re' : re.compile('>Query Chain:.+_'), 
        'templ_prot_aln_re' : re.compile('Query:.+'), 
        'query_prot_aln_re' : re.compile('Template:.+'), 
        }
    
    SMAP_dict = {}
    with open(SMAP_output, 'r') as output:
        for line in output:
            for regex in re_dict:
                if re_dict[regex].match(line):
                    if regex == 'templ_prot_re': 
                        templ_prot = re_dict[regex].findall(line)[0][-5:-1]
                        SMAP_dict[templ_prot] = None
                    elif regex == 'query_prot_re': 
                        query_prot = re_dict[regex].findall(line)[0][-5:-1]
                        SMAP_dict[query_prot] = None
                    elif regex == 'templ_prot_aln_re':
                        templ_prot_aln = re_dict[regex].findall(line)[0][10:] 
                        SMAP_dict[templ_prot] = templ_prot_aln.split(' ')
                    elif regex == 'query_prot_aln_re':
                        query_prot_aln = re_dict[regex].findall(line)[0][10:] 
                        SMAP_dict[query_prot] = query_prot_aln.split(' ')
    return SMAP_dict

def records_to_list(filename):
    from Bio import SeqIO
    records = []
    with open(filename, 'rU') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            records.append(record)
    return records

def write_data_to_csv(data, filename):
    import csv
    with open(filename, 'ab') as output:
        csv_obj = csv.writer(output, delimiter=',')
        csv_obj.writerows(data)

def read_csv(filename):
    import csv
    csvfile = open(filename, 'rb')
    csv_obj = csv.reader(csvfile, delimiter=',')
    return csv_obj

def dict_to_csv(dictfile, csv_name):
    from pandas import DataFrame

    seq_index_df = DataFrame.from_dict(dictfile)
    #print seq_index_df
    seq_index_df.to_csv(csv_name)


def clean_records(records_list):
    pass

def aln_merger(seq_1, seq_2):
    ''' Given two alignments for the same sequence, merges the two alignments
    and returns a single alignment '''
    from Bio.Seq import MutableSeq 

    mut_seq_1 = seq_1.tomutable()
    mut_seq_2 = seq_2.tomutable()

    for i_1, char_1 in enumerate(mut_seq_1):
        for i_2, char_2 in enumerate(mut_seq_2):
            if i_1 == i_2:
                if char_1 != char_2:
                    if char_1 == '-':
                        mut_seq_2.insert(i_2, '-')
                    elif char_2 == '-':
                        mut_seq_1.insert(i_1, '-')

    if mut_seq_1 == mut_seq_2:
        return mut_seq_1.toseq()
    elif len(mut_seq_1) > len(mut_seq_2):
        for num in range(len(mut_seq_1) - len(mut_seq_2)):
            mut_seq_1.pop()

        return mut_seq_2.toseq()
    elif len(mut_seq_2) > len(mut_seq_1):
        for num in range(len(mut_seq_2) - len(mut_seq_1)):
            mut_seq_2.pop()

        return mut_seq_2.toseq()

#seq1 = Seq('-AAA-G--F--')
#seq2 = Seq('ACAG-F-')
#
#merge_seq = aln_merger(seq1, seq2)
#print merge_seq
