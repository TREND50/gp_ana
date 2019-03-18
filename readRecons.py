# Script to plot teh results of reconstruction 
# Using files R[RunID]_planerecons.txt & R[RunID]_sphrecons.txt
# Produced wit the TREND recons software
# OMH March 18 2019

import os
import sys
import numpy as np
import pylab as pl
import plotAntennas

# Plane reconstruction analysis
pfile = "R"+sys.argv[1]+"_planerecons.txt"
a = np.loadtxt(pfile)
tux = a[:,1]
tux = tux-min(tux)
mult = a[:,2]
mults = np.unique(mult)
th = a[:,3]
phi = a[:,5]
th[th>90] = 180-th[th>90]

pl.figure(1)
pl.subplot(221)
pl.title('R{0} - Plane wave reconstruction'.format(sys.argv[1]))
for l in mults:
  pl.hist(th[mult==l],90,label='Mult={0}'.format(int(l)))
pl.xlim([0,90])
pl.xlabel('Zenith (deg)')
pl.legend(loc='best')
pl.subplot(222)
for l in mults:
  pl.hist(phi[mult==l],360,label='Mult={0}'.format(int(l)))
pl.xlim([0,360])
pl.xlabel('Azimuth (deg)')
pl.legend(loc='best')
pl.subplot(223)
for l in mults:
  pl.plot(tux[mult==l],th[mult==l],'+',label='Mult={0}'.format(int(l)))
pl.xlim([0,max(tux)])
pl.ylim([0,90])
pl.xlabel('Time (s)')
pl.ylabel('Zenith (deg)')
pl.legend(loc='best')
pl.subplot(224)
for l in mults:
  pl.plot(tux[mult==l],phi[mult==l],'+',label='Mult={0}'.format(int(l)))
pl.xlim([0,max(tux)])
pl.ylim([0,360])
pl.xlabel('Time (s)')
pl.ylabel('Azimuth (deg)')
pl.legend(loc='best')


# Spherical reconstruction analysis
sfile = "R"+sys.argv[1]+"_sphrecons.txt"
a = np.loadtxt(sfile)
mult = a[:,2]
xs = a[:,3]
ys = a[:,4]
zs = a[:,5]

# Select sources close to array only
close = (np.abs(xs)<2000) & (np.abs(ys)<4000)
xsc = xs[close]
ysc = ys[close]
zsc = zs[close]
mult = mult[close]
mults = np.unique(mult)

pl.figure(3)
pl.subplot(311)
pl.title('R{0} - Reconstructed source position'.format(sys.argv[1]))
for l in mults:
  pl.hist(xsc[mult==l],100,label='Mult={0}'.format(int(l)))
pl.legend(loc='best')
pl.xlabel('Easting (m)')
pl.subplot(312)
for l in mults:
  pl.hist(ysc[mult==l],100,label='Mult={0}'.format(int(l)))
pl.legend(loc='best')  
pl.xlabel('Northing (m)')
pl.subplot(313)
for l in mults:
  pl.hist(zsc[mult==l],100,label='Mult={0}'.format(int(l)))
pl.legend(loc='best') 
pl.xlabel('Altitude asl (m)')


pl.figure(2)
for l in mults:
  pl.plot(xsc[mult==l],ysc[mult==l],'+',label='Mult={0}'.format(int(l)))
pl.legend(loc='best') 
# Add antenna layout to plot
plotArray()

pl.show()

