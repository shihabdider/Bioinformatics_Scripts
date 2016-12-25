with open('pairwise_aa_substitution_vectors_self.csv') as csv:
    for line in csv:
        total = 0
        for num in line.split(','):
            total += float(num)
        print total
        
