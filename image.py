#!/usr/local/bin/python3

from matplotlib import pyplot as plt
import csv
data = []
with open('sout.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=';')
    for row in csv_reader:
        sub_data = []
        for elem in row:
            if elem != '':
                sub_data.append(float(elem))
        data.append(sub_data)
plt.imshow(data, interpolation='nearest')
plt.show()