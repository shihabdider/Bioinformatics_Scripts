from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Load records in fasta to a list
filename = 'top_ranked_residues.fasta'
def records_to_list(filename):
    records = []
    with open(filename, 'rU') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            records.append(record)
    return records

records = records_to_list(filename)

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

jaccard_set = list(read_csv('Jarccard_Similarity_Scores.csv'))

training_prots = []
for elem in jaccard_set:
    if elem[0] not in training_prots:
        training_prots.append(elem[0])
    if elem[1] not in training_prots:
        training_prots.append(elem[1])

print len(training_prots)

new_records = []
check_list = []
for prot in training_prots:
    for record in records:
        if prot in record.id:
            if prot not in check_list:
                check_list.append(prot)
                new_records.append(record)

with open('top_ranked_residues_training_set.fasta', 'w') as handle:
    for record in new_records:
        SeqIO.write(record, handle, 'fasta')
