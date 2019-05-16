# Script to plot the results of reconstruction
# Using files R[RunID]_planerecons.txt & R[RunID]_sphrecons.txt
# Produced wit the gp_recons software
# OMH March 2019

import os
import sys
import numpy as np
import pylab as pl
from plotAntennas import plotArray

c0 = 299792458
sigma_t = 15  #ns # Timing error
pl.ion()
pl.show()

DISPLAY = 1  # Switch to 1 for additionnal plots

#txtdir = "./testR222_without03/"
#txtdir = "./testR222/"
txtdir = "./"

def plot_recons(runid):

  # Load results of plane wave reconstruction
  pfile = txtdir+"R"+runid+"_planerecons_full.txt"
  if os.path.isfile(pfile) == False:
    print("No file",pfile,"! Generating it...") # Needed to compute Chi2 of recons
    loop_plot_delays(runid)

  a = np.loadtxt(pfile)
  tux = a[:,1]
  tux = tux-min(tux)
  mult = a[:,2]
  mults = np.unique(mult)
  th = a[:,3]
  phi = a[:,5]
  th[th>90] = 180-th[th>90]  # Mirror events reconstructed below ground above it
  chi2p = a[:,9]

  goodChi2p = (chi2p<1000)  # Select valid reconstruction only, set at CHi2<1000
  print(np.sum(goodChi2p),"valid plane recons out of",len(chi2p))
  th = th[goodChi2p]
  phi = phi[goodChi2p]
  tux = tux[goodChi2p]
  mult = mult[goodChi2p]
  mults = np.unique(mult)

  # Plot
  pl.figure(1)
  pl.subplot(221)
  pl.suptitle('R{0} - Plane wave reconstruction'.format(sys.argv[1]))
  for l in mults:
    pl.hist(th[mult==l],180,label='Mult={0}'.format(int(l)))
  pl.xlim([0,90])
  pl.xlabel('Zenith (deg)')
  pl.legend(loc='best')
  pl.subplot(222)
  for l in mults:
    pl.hist(phi[mult==l],720,label='Mult={0}'.format(int(l)))
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

  # Load results of spherical reconstruction analysis
  sfile = txtdir+"R"+sys.argv[1]+"_sphrecons_full.txt"
  if os.path.isfile(sfile) == False:
    print("No file",sfile,"! Generating it...")
    loop_plot_delays(runid)
  a = np.loadtxt(sfile)
  mult = a[:,2]
  xs = a[:,3]
  ys = a[:,4]
  zs = a[:,5]
  chi2s = a[:,9]

  # Select sources with valid recons & close to array only
  goodChi2s = (chi2s<1000)
  print(np.sum(goodChi2s),"valid spherical recons out of",len(chi2s))
  close = (np.abs(xs)<1e4) & (np.abs(ys)<1e4)
  xsc = xs[close & goodChi2s]
  ysc = ys[close & goodChi2s]
  zsc = zs[close & goodChi2s]
  mult = mult[close & goodChi2s]
  print(np.sum([close & goodChi2s]),"close and valid spherical recons out of",len(chi2s))
  mults = np.unique(mult)

  # Plot
  pl.figure(3)
  pl.subplot(311)
  pl.suptitle('R{0} - Reconstructed source position'.format(sys.argv[1]))
  for l in mults:
    pl.hist(xsc[mult==l],400,label='Mult={0}'.format(int(l)))
  pl.legend(loc='best')
  pl.xlabel('Easting (m)')
  pl.subplot(312)
  for l in mults:
    pl.hist(ysc[mult==l],400,label='Mult={0}'.format(int(l)))
  pl.legend(loc='best')
  pl.xlabel('Northing (m)')
  pl.subplot(313)
  for l in mults:
    pl.hist(zsc[mult==l],400,label='Mult={0}'.format(int(l)))
  pl.legend(loc='best')
  pl.xlabel('Altitude asl (m)')

  pl.figure(35)
  for l in mults:
    pl.plot(xsc[mult==l],ysc[mult==l],'+',label='Mult={0}'.format(int(l)))
  pl.legend(loc='best')

  # Add antenna layout to plot
  plotArray()

  pl.figure(4)
  pl.subplot(211)
  pl.hist(np.log10(chi2p+1),200)
  pl.xlabel('log10($\chi^2$ plane)')
  pl.subplot(212)
  pl.hist(np.log10(chi2s+1),200)
  pl.xlabel('log10($\chi^2$ sph)')
  # TBD: plot of Chi2 as a function of unit ID in event, in order to spot antenna with bad timing or bad position

def loop_plot_delays(runid):
  # Loop on all reconstructed coinc and call plot_delays in order to compute Chi2

  # First load ALL txt infos
  cfile = txtdir+"R"+runid+"_coinctable.txt"
  c = np.loadtxt(cfile)
  ccoincids = c[:,3]

  sfile = txtdir+"R"+sys.argv[1]+"_sphrecons.txt"
  s = np.loadtxt(sfile)
  scoincids = s[:,0]

  pfile = txtdir+"R"+sys.argv[1]+"_planerecons.txt"
  p = np.loadtxt(pfile)
  pcoincids = p[:,0]
  chi2ini = p[:,8]
  sel = (np.isinf(chi2ini)) | (np.isnan(chi2ini))
  p[sel,8] = 0

  uid, x, y, z =  getPos()

  # Then loop on coincs and build list of parameters for plot_delays()
  chi2sv, chi2pv = [],[]
  for i in scoincids:
    if int(i/100) == i/100:
      print("Processing recons",int(i),"in run",runid)
    ind = np.argwhere(ccoincids==i)
    ants = np.array(c[ind,1]).T[0]
    exp_delays = np.array(c[ind,4]).T[0]
    #
    ant_pos = []
    for a in ants:
      ind = np.argwhere(uid == a)[0] # Makes me crazy
      ant_pos.append([a, x[ind][0],y[ind][0], z[ind][0]])  #Northing, Westing, # Up
    ant_pos = np.array(ant_pos)
    #
    ind = np.argwhere(scoincids == i)[0]
    rec_source = np.array([s[ind,4],-s[ind,3],s[ind,5]]).T[0]
    #
    ind = np.argwhere(pcoincids == i)[0]
    rec_dir = [p[ind,3][0], p[ind,5][0]]

    # Now call plot_delays() and store resulting Chi2 value to array
    chi2s, chi2p = plot_delays(runid, int(i), exp_delays, rec_source = rec_source, rec_dir = rec_dir, ant_pos = ant_pos)
    chi2sv.append(chi2s)
    chi2pv.append(chi2p)

  # Append Chi2 value to recons results and write to file
  s = np.array(s)
  chi2sv = np.array(chi2sv)
  s = np.c_[s,chi2sv]  # Append chi2 column to result file
  sfile = "R"+sys.argv[1]+"_sphrecons_full.txt"
  np.savetxt(sfile,s,fmt='%f')  # Write to file
  p = np.c_[p,chi2pv]  # Append chi2 column to result file
  pfile = "R"+sys.argv[1]+"_planerecons_full.txt"
  np.savetxt(pfile,p,fmt='%f')  # Write to file


def plot_delays(runid, coincid, exp_delays = None, rec_source = None, rec_dir = None,ant_pos = None):
  # Plot expected vs measured trig delays

  # Load arrays of delays if not in argument
  if exp_delays is None:
    print("Processing recons",coincid," in run",runid)
    coincid = int(coincid)
    # Exp delays
    cfile = txtdir+"R"+runid+"_coinctable.txt"
    c = np.loadtxt(cfile)
    coincids = c[:,3]
    ind = np.argwhere(coincids==coincid)
    ants = np.array(c[ind,1]).T[0]
    exp_delays = np.array(c[ind,4]).T[0]
  else:
    ants = ant_pos[:,0]

  # Load antenna positions if not in argument
  if ant_pos is None:
    uid, x, y ,z = getPos()
    for a in ants:
      ind = np.argwhere(uid == a)[0] # Makes me crazy
      ant_pos.append([a, x[ind][0],y[ind][0], z[ind][0]])  #Northing, Westing, # Up
    ant_pos = np.array(ant_pos)

  # Load sph recons if not in argument
  if rec_source is None:
    sfile = txtdir+"R"+str(runid)+"_sphrecons.txt"
    a = np.loadtxt(sfile)
    coincids = a[:,0]
    ind = np.argwhere(coincids == coincid)[0]
    rec_source = np.array([a[ind,4],-a[ind,3],a[ind,5]]).T[0]

  #Now build expected delays for spherical hypothesis
  disv = ant_pos[:,1:4]-rec_source
  recs_delays = np.linalg.norm(disv,axis=1)/c0*1e9
  recs_delays = [recs_delays-np.min(recs_delays)][0]
  # Compute Chi2
  chi2s = np.sum(np.square(recs_delays-exp_delays)/(sigma_t*sigma_t))
  ndf = len(ants)-3
  chi2sndf = chi2s/ndf

  # Load plane recons if not in argument
  if rec_dir is None:
    pfile = txtdir+"R"+str(runid)+"_planerecons.txt"
    a = np.loadtxt(pfile)
    coincids = a[:,0]
    ind = np.argwhere(coincids == coincid)[0]
    rec_dir = [a[ind,3][0], a[ind,5][0]]

  # Now build expected delays for plane hypothesis
  # First build wavefront vector
  st = np.sin(np.radians(rec_dir[0]))
  ct = np.cos(np.radians(rec_dir[0]))
  sp =  np.sin(np.radians(rec_dir[1]))
  cp = np.cos(np.radians(rec_dir[1]))
  k = np.array([st*cp, st*sp,ct])
  # Now build delay vector
  disv = np.matmul(ant_pos[:,1:4],k)
  recp_delays = disv/c0*1e9
  recp_delays = recp_delays-np.min(recp_delays)
  # Compute Chi2
  chi2p = np.sum(np.square(recp_delays-exp_delays)/(sigma_t*sigma_t))
  ndf = len(ants)-2
  chi2pndf = chi2p/ndf

  # Display
  if DISPLAY and len(ants)>3:
  # and len(ants)>4 and chi2sndf>100:
  # and exp_delays[0]==103:
  # and len(ants)>4 and chi2sndf>100:
    print("Coinc",coincid,"R",runid,":")
    print("Experimental delays:",exp_delays)
    print("Reconstructed delays:",recs_delays)
    print("for point source reconstructed at location:",rec_source)
    print("Chi2/ndf =", chi2sndf)
    print("Reconstructed delays:",recp_delays)
    print("for reconstructed wave propag direction:",rec_dir)
    print("Chi2/ndf =", chi2pndf)

    pl.figure(1)
    pl.subplot(1,2,1)
    pl.plot(exp_delays,recp_delays,'sk')
    pl.title("Plane recons")
    pl.grid(True)
    pl.xlabel("Experimental trigger delays (ns)")
    pl.ylabel("Reconstructed trigger delays (ns)")
    xl = [-100, max(exp_delays),max(exp_delays)]
    pl.plot(xl,xl,'r')
    for i in range(len(ants)):
      pl.text(exp_delays[i]+100,recp_delays[i],str(int(ants[i])))
    textmax = max(max(exp_delays),max(recp_delays))
    pl.text(0,textmax*0.95,"($\Theta$,$\phi$)=({0},{1}) deg".format(round(rec_dir[0],1),round(rec_dir[1],1)))
    pl.text(0,textmax*0.9,"$\chi^2$/ndf = {0}".format(round(chi2pndf),1))

    pl.subplot(1,2,2)
    pl.plot(exp_delays,recs_delays,'sk')
    pl.plot(xl,xl,'r')
    pl.grid(True)
    pl.xlabel("Experimental trigger delays (ns)")
    pl.ylabel("Reconstructed trigger delays (ns)")
    for i in range(len(ants)):
      pl.text(exp_delays[i]+100,recs_delays[i],str(int(ants[i])))
    textmax = max(max(exp_delays),max(recs_delays))
    pl.text(0,textmax*0.95,"Source position=({0},{1},{2})m".format(round(rec_source[0],1),round(rec_source[1],1),round(rec_source[2],2)))
    pl.text(0,textmax*0.9,"$\chi^2$/ndf = {0}".format(round(chi2sndf,1)))
    pl.title("Spherical recons")
    pl.suptitle('Coinc {0} R{1}'.format(coincid,runid))

    input()
    pl.close('all')
  return chi2sndf, chi2pndf

if __name__ == '__main__':
  loop_plot_delays(sys.argv[1])
  if len(sys.argv)==3:
    plot_delays(sys.argv[1],sys.argv[2])

  plot_recons(sys.argv[1])
