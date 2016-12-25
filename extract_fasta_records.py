from Bio import SeqIO, Seq

filename = 'gross-alignment.fasta'
def records_to_list(filename):
    records = []
    with open(filename, 'rU') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            records.append(record)
    return records

records = records_to_list(filename)

for record in records:
    filename = record.id
    output = "".join([c for c in filename if c.isalpha() or c.isdigit() or
        c==' ']).rstrip()
    with open('bs_{0}.fasta'.format(output), 'w') as handle:
        SeqIO.write(record, handle, 'fasta')
