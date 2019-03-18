# Plot GP35 antenna array from ants.txt
# OMH January 2019

import numpy as np
import pylab as pl


def plotArray():
  # Load antenna position file
  antfile = "ants.txt"
  pos = np.loadtxt(antfile,delimiter=',')
  uid=pos[:,0]
  x=pos[:,2]
  y=pos[:,1]
  z=pos[:,3]

  # Status of antennas
  off=[16,22,19]
  dead=[25]
  notsure=[3,6]
  isin = np.setxor1d(uid,off)
  isin = np.setxor1d(isin,dead)
  isin = np.setxor1d(isin,notsure)

  # Plot
  pl.figure(2)
  pl.plot(x[np.in1d(uid,off)],y[np.in1d(uid,off)],'ok')
  pl.plot(x[np.in1d(uid,dead)],y[np.in1d(uid,dead)],'xr')
  pl.plot(x[np.in1d(uid,notsure)],y[np.in1d(uid,notsure)],'vb')
  pl.plot(x[np.in1d(uid,isin)],y[np.in1d(uid,isin)],'og')

  for i in range(len(uid)):
    pl.text(x[i]+10,y[i],str(int(uid[i])))

  pl.grid(True)
  pl.xlabel("Easting (m)")
  pl.ylabel("Northing (m)")
  pl.gca().set_aspect('equal')
  pl.show()


if __name__ == '__main__': 
  plotArray();
