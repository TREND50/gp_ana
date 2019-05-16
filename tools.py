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

def getPos():
# Reads GP35 antenna positions in ants.csv, returns coordinates and dumps them to file, using GRAND conventions
  antfile = "ants.csv"
  data = getCSVData(antfile)
  uid = np.array(data["unit ID"])
  x = data["northing"]  # x: Northing
  x = np.array(list(map(float, x)))
  y = data["easting"]  # y: Westing
  y = -np.array(list(map(float, y)))
  z = data["elevation"]
  z = np.array(list(map(float, z)))

  # Select positions where units are deployed
  indb = []
  for i, j in enumerate(uid):
     if 'b' in j:  # Letter 'b' has to appear in field "unit ID" of CSV file
         indb.append(i)
  uid = uid[indb]
  x = x[indb]
  y = y[indb]
  z = z[indb]

  uid = list(map(lambda each:each.strip("b"), uid))
  uid = np.array(list(map(int, uid)))

  res = np.hstack((uid, x, y , z))
  res = np.reshape(res,(4,len(uid)))
  res = res.T  # Insn't this crazy???
  np.savetxt("ants.tmp", res, fmt='%d %.3f %.3f %.3f')
  return [uid, x, y , z]

def dumpCoord_recons():
# Reads GP35 antenna positions in ants.csv, dumps them to file for reconstruction (ie using TREND conventions)
# Fetch coordinates
  [uid, x, y , z] = getPos()

  res = np.zeros(shape=(35,4))
  res[:,0] = range(1,36)  # Unit ID
  res[uid-1,1] = -y # Easting
  res[uid-1,2] = x # Northing
  res[uid-1,3] = z

  print(res)
  np.savetxt("positions_GP35.tmp", res, fmt='%d %.3f %.3f %.3f')
