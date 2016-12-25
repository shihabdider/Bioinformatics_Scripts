'''Extracts specific columns based on an index (either a gapped index, or
non-gapped index)'''

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Load records in fasta to a list
#filename = 'gross-alignment.fasta'
filename = 'HRH1_EDNRB_FATCAT.fasta'
def records_to_list(filename):
    records = []
    with open(filename, 'rU') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            records.append(record)
    return records

records = records_to_list(filename)
training_prots_records = records_to_list(
    'top_ranked_residues_training_set.fasta')
# Extract index of columns 
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

#sorted_features = list(read_csv('sorted_features.csv'))
#absolute_index = []
#for num in range(10):
#    absolute_index.append(int(sorted_features[num][2]))
#
#print absolute_index
absolute_index = []

bs_residue_index = [6346, 6612, 5991, 6085, 6080, 6682, 6090, 6132,
        5920, 6333]
binding_site_5HT1B = [129, 130, 133, 134, 199, 200, 201, 216, 327,
        330, 334, 348, 351, 352, 355, 359] 
binding_site_4NC3 = [62, 353, 395, 396, 397]

binding_site_4MBS = [37, 86, 89, 108, 109, 182, 195, 248, 251, 259, 283]
binding_site_HRH1 = [107, 108, 111, 112, 158, 428, 431, 432, 458]
binding_site_adjust_HRH1 = [pos-27 for pos in binding_site_HRH1]
#with open('agonist2_site_2YDO.txt', 'r') as handle:
#    for line in handle:
#         bs_residue_index.append(int(line.split()[2]))

def extract_column_index(records, residue_index, get_columns=True,
        prot_template_id=None):
    global training_prots_records
    seq_index_dict = {}
    if get_columns:
        for record in records:
            if prot_template_id in record.id:
                print record.id, 'Extracting specified residues...'
                total_count = -1
                seq_count = 1
                
                site_list = []
                for char in record.seq:
                    total_count += 1
                    if char != '-':
                        if seq_count in residue_index:
                            print seq_count
                            print total_count
                            #absolute_index.append(total_count)
                            print record.seq[total_count]
                            site_list.append([total_count, seq_count, char])
                        seq_count += 1

                seq_index_dict[record.id] = site_list

    else:
        #for train_record in training_prots_records:
            for record in records:
                #if train_record.id == record.id:
                print record.id, 'Extracting specified residues...'

                total_count = -1
                seq_count = 1
                
                site_list = []
                for char in record.seq:
                    total_count += 1
                    if char != '-':
                        if total_count in residue_index:
                            print seq_count
                            print total_count
                            #absolute_index.append(total_count)
                            print record.seq[total_count]
                            site_list.append([total_count, seq_count, char])
                    seq_count += 1
                    
                seq_index_dict[record.id] = site_list

    return seq_index_dict

seq_index_dict = extract_column_index(records, binding_site_adjust_HRH1,
        prot_template_id='hrh1')

#print seq_index_dict

key_id = []
for key in seq_index_dict:
    if 'hrh1' in key:
        residue_index = seq_index_dict[key]

gapped_index = [pos[0] for pos in residue_index]    
print gapped_index

seq_gapped_index_dict = extract_column_index(records, gapped_index,
        get_columns = False)

print seq_gapped_index_dict
def dict_to_csv(dictfile, csv_name):
    from pandas import DataFrame

    seq_index_df = DataFrame.from_dict(dictfile)
    #print seq_index_df
    seq_index_df.to_csv(csv_name)

#dict_to_csv(seq_index_dict, 'HRH1_EDNRB_binding_site_index.csv')

## Extract columns using the index
#
#new_records = []
#for record in records:
#    bs_residues = ''
#    for i in absolute_index:
#        bs_residues += record.seq[i]
#    new_record = SeqRecord(Seq(bs_residues), id=record.id)
#    new_records.append(new_record)
#
## Write residues to file
#with open('top_ranked_residues.fasta', 'w') as handle:
#    for record in new_records:
#        SeqIO.write(record, handle, 'fasta')
