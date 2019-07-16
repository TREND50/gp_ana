# Plot GP35 antenna array from ants.csv
# OMH May 2019

import numpy as np
import pylab as pl
from tools import getCSVData

def plotArray():
  # Load antenna position file
  antfile = "ants.csv"
  data = getCSVData(antfile)
  uid = np.array(data["unit ID"])
  status = np.array(data["status"])
  x = data["easting"]  # Westing ==> x
  x = np.array(list(map(float, x)))
  y = data["northing"]  # Northing ==> y
  y = np.array(list(map(float, y)))
  z = data["elevation"]
  z = np.array(list(map(float, z)))

  # Select positions where units are deployed
  indb = []
  for i, j in enumerate(uid):
     if 'b' in j:  # Letter 'b' has to appear in field "unit ID" of CSV file
         indb.append(i)
  uid = uid[indb]
  status = status[indb]
  x = x[indb]
  y = y[indb]
  z = z[indb]

  # Status
  ok = [status == "ok"]
  nocom = [status == "no ping"]
  dead = [status == "no signal"]
  part = [status == "x"]
  bad = [status == "noisy"]
  cal = [status == "calibrator"]

  # Plot
  pl.figure(35)
  pl.plot(x[tuple(nocom)],y[tuple(nocom)],'.k',label = "no com")
  pl.plot(x[tuple(dead)],y[tuple(dead)],'xr',label="no signal")
  pl.plot(x[tuple(part)],y[tuple(part)],'ob',label="partial")
  pl.plot(x[tuple(bad)],y[tuple(bad)],'or',label="noisy")
  pl.plot(x[tuple(ok)],y[tuple(ok)],'og',label="ok")
  pl.plot(x[tuple(cal)],y[tuple(cal)],'hm',label="emitter")

  for i in range(len(uid)):
    pl.text(x[i]+10,y[i],uid[i])

  pl.grid(True)
  pl.xlim([-500,+500])
  pl.xlabel("Easting (m)")
  pl.ylabel("Northing (m)")
  pl.legend(loc="best")
  pl.title('GP35 - June 30, 2019')
  pl.gca().set_aspect('equal')
  pl.show()


if __name__ == '__main__':
  plotArray();
