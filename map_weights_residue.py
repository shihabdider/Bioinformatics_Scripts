import csv

def write_data_to_csv(data, filename):
    with open(filename, 'ab') as output:
        csv_obj = csv.writer(output, delimiter=',')
        csv_obj.writerows(data)

def read_csv(filename):
    csvfile = open(filename, 'rb')
    csv_obj = csv.reader(csvfile, delimiter=',')
    return csv_obj

feature_weights = []
cleaned_pos_index = []
for item in read_csv('cleaned_pos_index.csv'):
    for elem in item:
        cleaned_pos_index.append(int(elem))

with open('feature_weights_new.txt', 'r') as weights:
    for i, line in enumerate(weights):
        feature_weights.append([float(line[:-1]), i, cleaned_pos_index[i]])


sorted_features = sorted(feature_weights, key = lambda x: float(x[0]),
        reverse=True)

print sorted_features

write_data_to_csv(sorted_features, 'sorted_features.csv')
