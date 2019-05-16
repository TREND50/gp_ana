import csv
import numpy as np
import sys
import pylab as pl
import os

def getCSVData(filename):
# Reac CSV file and return result as dict
    vals = []
    print('Reading',filename)
    with open(filename) as csvfile:
        reader = csv.reader(csvfile)
        ind = 0
        for row in reader:
            if ind == 0:  # Field definition
                dim = len(row)
                keys = row
                print("There are",dim,"fields in file",filename)
                print(keys)
            else:
                vals.append(row)
            ind +=1

# Now reorganize data in dict
    vals = np.array(vals)
    values = [vals[:,col].tolist() for col in range(dim)]  # Data as list of columns to fit in dict
    data = dict(zip(keys, values))
#    data = dict.fromkeys(keys, values)  # Does not work
#    print(data["latitude"])
    return(data)
