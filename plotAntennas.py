import numpy as np
import pylab as pl

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
isin = np.setxor1d(uid,off)
isin = np.setxor1d(isin,dead)

# Plot
pl.figure(1)
pl.plot(x[np.in1d(uid,off)],y[np.in1d(uid,off)],'ok')
pl.plot(x[np.in1d(uid,dead)],y[np.in1d(uid,dead)],'xr')
pl.plot(x[np.in1d(uid,isin)],y[np.in1d(uid,isin)],'og')
for i in range(len(uid)):
  pl.text(x[i]+10,y[i],str(int(uid[i])))

pl.grid(True)
pl.xlabel("Easting (m)")
pl.ylabel("Northing (m)")
pl.show()