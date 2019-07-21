import os
import sys
import numpy as np
import pylab as pl

import readData as rd


if len(sys.argv)!=3:
    print("Usage: >loopRead RUNID RUNTYPE")
    sys.exit()

nrun = sys.argv[1]
typ = sys.argv[2]


#First load run
pyf = rd.load_data(nrun,typ)

# Now build stat
#muba, sigba, muaa = np.zeros(shape=(len(antsIn),4)),np.zeros(shape=(len(antsIn),4)),np.zeros(shape=(len(antsIn),4))
ants = range(1,36)
antsIn = []
muba,sigba,muaa = [],[],[]
for i,ant in enumerate(ants):
    try:
        m,s,a = rd.display_events(pyf=pyf,typ=typ,tid=ant)
        muba.append(m[0])
        sigba.append(s[0])
        muaa.append(a[0])
        antsIn.append(ant)
    except SystemExit:
        print("No data for unit b",ant)

# Now display
muba = np.array(muba)
sigba = np.array(sigba)
muaa = np.array(muaa)
antsIn = np.array(antsIn)
lab = ["X","Y","Z"]
pl.figure(1)
for k in range(3):
    pl.subplot(131)
    pl.plot(antsIn,muaa[:,k],'o',label=lab[k])
    pl.xlabel("Antenna ID")
    pl.ylabel("$\mu_{Amp}$ (V$_{ADC}$)")
    pl.grid(True)
    pl.subplot(132)
    pl.plot(antsIn,muba[:,k],'o',label=lab[k])
    pl.xlabel("Antenna ID")
    pl.ylabel("$\mu_{bline}$ (V$_{ADC}$)")
    pl.legend(loc="best")
    pl.grid(True)
    pl.subplot(133)
    pl.plot(antsIn,sigba[:,k],'o',label=lab[k])
    pl.xlabel("Antenna ID")
    pl.ylabel("$\sigma_{bline}$ (V$_{ADC}$)")
    pl.grid(True)
pl.suptitle('{0}{1}'.format(typ,nrun))
pl.show()
