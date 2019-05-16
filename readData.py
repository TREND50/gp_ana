# Script to read GP35 data produced by the RUST DAQ software
# and in particular build the coinctable.txt file
# to be used with the gp_recons software to reconstruct shower direction of origin

# OMH January 2019



import os
import sys
sys.path.append("../")
import pyef

import numpy as np
import pylab as pl
import yaml
import scipy
from scipy.optimize import curve_fit

DISPLAY = 0  # Switch to 1 for additionnal plots
ULASTAI = 0
if ULASTAI:
  datafolder = "/mnt/disk"
else:  # @ local computer
  datafolder = "/home/martineau/GRAND/GRANDproto35/data/ulastai"

pl.ion()
c0 = 299792458
utcSLC = []
maxCoarse = []

def twos_comp(val, bits):
    """compute the 2's compliment of int value val"""
    if (val & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
        val = val - (1 << bits)        # compute negative value
    return val


def loadMaxCoarse(runid):
  # Load max coarse info from SLC file
  global IDsin
  allIDs = []
  allUTC = []
  allMaxCoarse = []

  # First check that file exists
  slcfile = datafolder+"/S"+runid+".yaml"
  if os.path.isfile(slcfile) is False:
     print('File ',slcfile,'does not exist. Aborting.')
     IDsin = []
     return

  # Now read SLC file and dump infos into numpy arrays
  print('Scanning SLC file',slcfile)
  print("Loading SLC data...")
  dataf=yaml.load_all(open(slcfile))
  print("Done.")
  for d in dataf:
    if d['msg_type']=='SLC':
      allIDs.append(d['source_ip'][3]-100)
      allUTC.append(d['received_timestamp'][0])
      allMaxCoarse.append(d['max_coarse'])

  allIDs = np.array(allIDs)
  allUTC = np.array(allUTC)
  allMaxCoarse = np.array(allMaxCoarse)
  IDsin = np.unique(allIDs)
  print("MaxCoarse info retrieved for following units:",IDsin)

  # Now order info according to Unit ID
  for uid in IDsin:
    ind = np.nonzero(allIDs==uid)
    utcSLC.append(allUTC[ind])
    maxCoarse.append(allMaxCoarse[ind])


def getMaxCoarse(uid,utcsec):
  # Retrieve maxCoarse value closests in time to event triggered @ t = utcsec
  i = np.nonzero(IDsin == uid)[0]
  if len(i)>0:  # Unit found in SLC data
    i = i[0]
    indt = np.argmin(np.abs(utcSLC[i]-utcsec))
    #print(utcsec,utcSLC[i][indt],maxCoarse[i][indt])
    # TBD: implement "closest time" condition (<1h)?
  else:  # Unit not found in LSC data
    return 0

  return maxCoarse[i][indt]


def build_distmat():
  # Build matrix of d(ant1,ant2)
  antfile = "ants.txt"
  pos = np.loadtxt(antfile,delimiter=',')
  a = np.array([('a', 2), ('c', 1)], dtype=[('x', 'S1'), ('y', int)])
  uid=pos[:,0]
  nants = len(uid)
  x=pos[:,2]
  y=pos[:,1]
  z=pos[:,3]
  #p = [x, y, z]  # Why cannot numpy build matrixes with this syntax?????????? Makes me crazy
  p = np.vstack((x,y,z))
  d = np.ndarray(shape=(nants,nants))
  for i in range(nants):
    for j in range(nants):
      #print(np.linalg.norm(p[:,j]-p[:,i]))
      d[i,j] = np.linalg.norm(p[:,j]-p[:,i])
      d[j,i] = d[i,j]

  return uid,d


def build_coincs(trigtable,uid,d):
  # Search for coincs
  nrun = sys.argv[1]
  print("Now searching for coincidences...")
  ntrigs = np.shape(trigtable)[0]
  tmax = np.max(d)/c0*1e9*1.1  # Factor 1.1 to give some flexibility
  # TBD: adjust tmax for each pair of antennas
  uid = trigtable[:,0]  # Vector of unit UDs
  secs = trigtable[:,1]  # Vector of seconds info
  secscor = secs-min(secs)  #  Use first second as reference. Otherwise "times" field is too long and subsequent operations fail...
  nsecs = trigtable[:,2] # Vector of nanoseconds info
  times = secscor*1e9+nsecs # Build complete time info. Units = ns.

  #
  i = 0
  coinc_nb = 0
  delays = []
  uid_delays = []
  filename = 'R{0}_coinctable.txt'.format(sys.argv[1])  # File where coincs should be written, latter used for source reconstruction
  while i<ntrigs:
  # Loop on all triggers in table
    trig_ref = times[i]
    id_ref = uid[i]
    tsearch = times[i:-1]-trig_ref
    tsearch = tsearch[np.argwhere(tsearch<tmax)].T[0]  #  Search in causal timewindow. Transpose needed to get a line vector and avoid []
    idsearch = uid[i:i+len(tsearch)]
    #print(i,"*** Reference unit:",i,id_ref,uid[i],secs[i],nsecs[i],trig_ref,times[i],", now looking for coincs within",tmax,"ns")
    _, coinc_ind = np.unique(idsearch, return_index = True) # Remove multiple triggers from a same antenna
    if len(coinc_ind)>3:  # Requires 4 antennas at least to perform recons
      # there are events in the causal timewindow
      coinc_nb = coinc_nb+1  # INcrement coinc counter
      #print(i,"*** Reference unit:",id_ref,times[i],", now looking for coincs within",tmax,"ns")
      #print(np.argwhere(tsearchini<tmax),tsearchini[0:10])
      #print("Units in causal timewindow:",i,i+len(tsearch),idsearch,tsearch,secs[i:i+len(tsearch)],nsecs[i:i+len(tsearch)])
      #print("From different units:",idsearch[others],tsearch[others],others)
      mult = len(coinc_ind)
      print(coinc_nb,": possible coinc at (",int(secs[i]),"s;",nsecs[i],"ns) between",mult,"units:",idsearch[coinc_ind],tsearch[coinc_ind])

      coinc_ids = idsearch[coinc_ind]
      coinc_ids = np.array([coinc_ids])  # Anybody able to explain why coinc_id is not a numpy.array???
      coinc_times = tsearch[coinc_ind]
      coinc_times = np.array([coinc_times])

      # Write to file
      # Format: Unix sec; Unit ID, Evt Nb, Coinc Nb, Trig time (ns)
      evts = np.array([range(i,i+mult)],dtype = int).T
      one = np.ones((mult,1),dtype=int)
      this_coinc = np.hstack((secs[i]*one,coinc_ids.T,evts,coinc_nb*one,coinc_times.T))
      if coinc_nb == 1:
        all_coincs = this_coinc
      else:
        all_coincs = np.concatenate((all_coincs,this_coinc))

      # Now load delay info (only for histos)
      delays = np.concatenate((delays,coinc_times[coinc_ids!=id_ref]))
      uid_delays = np.concatenate((uid_delays,coinc_ids[coinc_ids!=id_ref]))
      i = i+len(coinc_ind)

    else:
      i = i+1

  np.savetxt(filename,all_coincs,fmt='%d')  # Write to file

  uid_delays = np.array(uid_delays)
  uid = np.unique(uid_delays)
  pl.figure(8)
  for i in uid:
    h,b,_ = pl.hist(delays[uid_delays==i],200,label='ID{0}'.format(int(i)))
  pl.xlim([0,tmax])
  pl.legend(loc='best')
  pl.xlabel('Trigger delay (ns)')
  pl.savefig('delays_R{0}'.format(sys.argv[1]))
  filename = "R{0}_trig_delays.npz".format(nrun)
  np.savez("R{0}_trig_delays".format(nrun),delays,uid_delays)  # Write delay histogram to file for faster access for further work on this distribution (see fitDelays())

  if DISPLAY:
    pl.show()


def fitDelays(filename):
  # Fit trigger delay distribution with gaussian distrib in order to estimate timing resolution
  a = np.load(filename)
  delays = a['arr_0']
  uids = a['arr_1']
  uid = np.unique(uids)

  pl.figure(12)
  for i in uid:
        h,b,_ = pl.hist(delays[uids==i],1000,label='ID{0}'.format(int(i)))

	# Perform Gaussian fit
        valMax = np.max(h)
        if valMax>10:
          posMax = b[np.argmax(h)]
          p0 = [valMax,posMax,20.]
          bc = (b[:-1] + b[1:])/2
          coeff, var_matrix = curve_fit(gauss, bc, h, p0=p0)
          h_fit = gauss(bc, *coeff)
          #pl.plot(bc, h, label='ID{0}'.format(int(i)))
          pl.plot(bc, h_fit)
          print('Unit',int(i),': mean trig delay',coeff[1],'ns latter than 1st trigger')
          print('Trigger time dispersion = ', coeff[2]/np.sqrt(2),'ns')

  pl.xlim([0,6000])
  pl.legend(loc='best')
  pl.xlabel('Trigger delay (ns)')
  input()

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


def get_time(nrun=None,pyf=None):
  # Retrieves time info from datafile and order it in increasing order
  print("Now building trigger time table from data...")
  if pyf == None:
    print("No pyef object, loading it from run number.")
    if nrun == None:
      print("get_time error! Pass run number or pyef object as argument")
      return
    pyf = load_data(nrun)

  if nrun == None:
    nrun = sys.argv[1]

  nevts = len(pyf.event_list)
  # Loop on all events
  secs = []
  nsecs = []
  ttimes = []
  IDs = []
  i = 0
  for evt in f.event_list:
    #print("\n\n!!!New event!!!")
    # Loop on all units in the event (at present should be only one)
    for ls in evt.local_station_list:
      # Loop on all units involved in event (should be one only at this stage)
      #print("Local unit info")
      #ls.display()
      uid = int(ls.header.ls_id-356) # 16 lowest bits of IP adress --256 to go down to 8 lowest digits of IP adress & -100 to go down to antenna ID
      #IDs.append(uid)
      #nsec = ls.header.gps_nanoseconds  # GPS info for that specific unit

    #if uid == 3: # Dirty trick to exclude one antenna from analysis. Redondant with event.header.event_nsec at this stage (one unit per event only)
    #  print("Skipping unit 03")
    #  continue

    if uid == 20:  # Dirty trick to exclude one antenna from analysis
        #  print("Skipping unit 03")
        continue
    IDs.append(uid)
    h = evt.header
    sec = h.event_sec
    nsec = h.event_nsec
    # Now correct time from maxCoarse value
    maxcoarse = getMaxCoarse(uid,sec)
    cor=125e6/(getMaxCoarse(uid,sec)+1)
    if abs(cor-1)>0.01: # Abnormal correction value
      cor = 1
    nsec = nsec*cor
    nsecs.append(nsec)
    secs.append(sec)

    print("Event",i,", ID=",uid,",Time=",sec,nsec)  # Can be used to check that SSS is OK.
    #print("Event info")
    #evt.display()
    #print("Header info")
    #h.display()
    i = i+1

  secs = np.array(secs)
  nsecs = np.array(nsecs)
  # Build total time info. Warning: set 1st event as reference otherwise value too large and argsort fails!!!
  if min(secs)<1:
    print("Error!!! Minimal second info =",min(secs),". Abort.")
    return
  # Build time info
  ttimes = (secs-min(secs))*1e9+nsecs  # Set 1st event as reference
  ttimes = np.array(ttimes,dtype=int)
  ttimes = (ttimes-min(ttimes))/1e9  # in seconds
  # Order in time
  IDs = np.array(IDs)
  units = np.unique(IDs)
  ind = np.argsort(ttimes)
  IDs_ordered = IDs[ind]
  secs_ordered = secs[ind]
  nsecs_ordered = nsecs[ind]
  res = np.vstack((IDs_ordered,secs_ordered,nsecs_ordered))
  res = res.T

  # Check for errors
  # TBD: add  flag for events with time info = 0 (ie no GPS info)

  # Check for time offsets
  #tdif = np.diff(secs)
  #aid = np.argwhere(tdif<0)
  #for i in aid:
  #   print("***Warning! Possible error on SSS value for board",IDs[i+1],":\nevent",i," on board",IDs[i],": SSS =",secs[i],"\nevent",i+1," on board",IDs[i+1],": SSS =",secs[i+1],"\nevent",i+2," on board",IDs[i+2],": SSS =",secs[i+2])

  # Check for time jumps in the past
  for uid in units:
      tdif = np.diff(secs[IDs==uid])
      aid = np.argwhere(tdif<0)
      for i in aid:
        #print(ind)
        #i = ind[1]
        #print(i,j,aid[j],xx)
        print("***Warning! Jump in past for unit",uid,"from SSS =",secs[IDs==uid][i],"to SSS =",secs[IDs==uid][i+1])

  # Plot a few graphs to check run quality
  pl.figure(1)
  pl.subplot(211)
  for uid in units:
    pl.plot(nsecs[IDs==uid],label=uid)
  pl.xlabel('Event nb')
  pl.ylabel('Nanosec counter value')
  pl.legend(loc='best')
  pl.subplot(212)
  for uid in units:
    pl.hist(nsecs[IDs==uid],100,label=uid)
  pl.xlabel('Nanosec counter value')
  pl.xlim([0,1e9])
  pl.legend(loc='best')

  pl.figure(2)
  for uid in units:
    pl.plot(ttimes[IDs==uid],label=uid)
  pl.xlabel('Event nb')
  pl.ylabel('Trigger time (s)')
  pl.legend(loc='best')
  pl.title('Event rate')
  pl.savefig('EventRate_R{0}'.format(sys.argv[1]))

  pl.figure(3)
  pl.hist(IDs,100)
  pl.xlabel("Unit ID")
  pl.ylabel("Nb of events")
  pl.savefig('NbEvents_R{0}'.format(sys.argv[1]))

  for uid in units:
    pl.figure(uid*100)
    delta_trig = np.diff(ttimes[IDs==uid])*1000 # in ms
    pl.hist(delta_trig[ (delta_trig<200) & (delta_trig>0)],1000)
    pl.title('Unit {0}'.format(uid))
    pl.xlabel("$\Delta$t trig (ms)")

  dur = max(ttimes)
  print(nevts,"events in run",nrun)
  print("Run duration:",dur,"seconds.")
  print("Units present in run:")
  for uid in units:
    print("Unit",uid,":",np.shape(np.where(IDs==uid))[1],"events.")
  input()

  if DISPLAY:
    pl.show()

  # Format: ID sec nsec (ordered by increasing time)
  return res



def display_events(nrun=None,pyf=None,typ="R",tid=None):
  # Display timetraces for unit tid
  print("Now assembling stats for timetraces of unit",tid)
  lab = ['X','Y','Z','Cal']
  if pyf == None:
    print("No pyef object, loading it from run number.")
    if nrun == None:
      print("display_events() error! Pass run number or pyef object as argument")
      return
    pyf = load_data(nrun,typ)

  offset = 90
  if typ == "C":
      nCh = 4
      win = range(3,2*offset) # Skip 1st 3 points because could be left overs from previous event
  else:
      nCh = 3
      win = range(3,offset-10) # Skip 1st 3 points because could be left overs from previous events

  if nrun == None:
    nrun = sys.argv[1]
  nevts = len(pyf.event_list)
  mub,sigb,imax,Amax = [],[],[],[]
  j = 0
  for evt in f.event_list:
    j = j+1
    if j/1000 == int(j/1000):
      print("Processing event",j,"/",nevts)

    # Loop on all events
    for ls in evt.local_station_list:
      # Loop on all units involved in event (should be one only at this stage of software)
      uid = ls.header.ls_id - 356 # Remove 255 for digits 8-16 in IP adress & 100 for unit ID
      if tid == None or uid == tid:  #  display events for unit tid (or all events from all units if tid == None)
        # Access data
        raw = ls.adc_buffer
        hraw = [hex(int(a)) for a in raw]  # Transfer back to hexadecimal
        draw = [twos_comp(int(a,16), 12) for a in hraw] #2s complements
        draw = np.array(draw)*1./2048  # in Volts
        nsamples = int(ls.header.trace_length/4)  # draw corresponds to 4 channels
        #offset = int(nsamples/2.0)  # Offset position at center of waveform
        #print nsamples,"samples per channel --> offset = ",offset
        thisEvent = np.reshape(draw,(4,nsamples));
        thisEvent = pow(10,(thisEvent+np.min(thisEvent)))  # Correct for logarithmic amplification of power detector. Then


        if DISPLAY:
            pl.figure(1)
            tmus = np.array(range(nsamples))*20e-3  # Time axis in µs
            evtnb = ls.header.event_nr
            for i in range(nCh):
                pl.plot(tmus[3:],thisEvent[i][3:],label=lab[i])
                # Draw line at expected trigger position
                #pl.plot([tmus[int(nsamples/2)+15], tmus[int(nsamples/2)+15]],[np.min(thisEvent[:][:]),np.max(thisEvent[:][:])])
                pl.title('R{0} Evt{1} Antenna {2}'.format(nrun,evtnb,uid))
                pl.xlim(tmus[3],max(tmus))
                pl.xlabel('Time ($\mu$s)')
                pl.ylabel('10$^{Voltage}$')
                pl.legend(loc="best")
            pl.show()
            input()
            pl.close('all')

        # Build stats
        imub,isigb,iimax,iAmax = np.zeros((nCh,1)),np.zeros((nCh,1)),np.zeros((nCh,1)),np.zeros((nCh,1))
        for i in range(nCh):
          imub[i] = np.mean(thisEvent[i][win])  # Restrict to baseline
          isigb[i] = np.std(thisEvent[i][win]) # Restrict to baseline
          iimax[i] = int(np.argmax(thisEvent[i][3:])+3)
          iAmax[i] = thisEvent[i][int(iimax[i])]

        #print(imub,isigb,iimax,iAmax)
        imax.append(iimax)
        Amax.append(iAmax)
        mub.append(imub)
        sigb.append(isigb)

  # Now plot full stat infos
  imax = np.array(imax)
  Amax = np.array(Amax)
  mub = np.array(mub)
  sigb = np.array(sigb)
  for k in range(nCh):
      if 1:
          good = np.sum( (imax[:,k]>104) & (imax[:,k]<108))
          abline = np.sum( (Amax[:,k]<0))
          azero = np.sum( (Amax[:,k]==0))
          print('Channel',k,': good events=',good,'/',np.shape(imax)[0],'=',float(good)/np.shape(imax)[0])
          print('Channel',k,': Max at zero=',azero,'/',np.shape(imax)[0],'=',float(azero)/np.shape(imax)[0])
          print('Channel',k,': Max < zero=',abline,'/',np.shape(imax)[0],'=',float(abline)/np.shape(imax)[0])

      pl.figure(tid*100+21+k)
      pl.subplot(231)
      pl.hist(mub[:,k],offset*2)
      pl.xlabel('Baseline mean')
      #pl.title('Board {0}'.format(tid))
      pl.grid(True)

      pl.subplot(235)
      pl.plot(mub[:,k],'+')
      pl.plot(Amax[:,k],'o')
      pl.xlabel('Event ID')
      pl.ylabel('Mean amp (bline & max)')
      #pl.title('Board {0}'.format(tid))
      pl.grid(True)

      pl.subplot(234)
      pl.xlabel('Index of signal max')
      pl.hist(imax[:,k],offset*2)
      #pl.title('Board {0}'.format(tid))
      pl.grid(True)

      pl.subplot(236)
      pl.xlabel('Max amplitude')
      pl.hist(Amax[:,k],offset*2)
      #pl.title('Board {0}'.format(tid))
      pl.grid(True)

      pl.subplot(232)
      diffAmp = Amax[:,k]-mub[:,k]
      pl.hist(sigb[:,k],offset*2)
      pl.xlabel('Bline std dev')
      #pl.title('Board {0}'.format(tid))
      pl.grid(True)

      print('Channel',k,': Peak @ ',np.mean((Amax[:,k][0])),'V, std dev=',np.std((Amax[:,k])),'V, rel error=',np.std((Amax[:,k]))/np.mean((Amax[:,k]))*100,'%')
      print('Channel',k,': Peak - bline @ ',np.mean((diffAmp)),'V, std dev=',np.std((diffAmp)),'V, rel error=',np.std((diffAmp))/np.mean((diffAmp))*100,'%')

      pl.subplot(233)
      pl.plot(mub[:,k],sigb[:,k],'+')
      pl.xlabel('Baseline mean')
      pl.ylabel('Bline std dev')
      #pl.title('Board {0}'.format(tid))
      pl.grid(True)

  print("Done. If no event was displayed,then this means that there was no event recorded for unit",tid,"in run",nrun,".")

def load_data(nrun,typ="R"):
  # Loads run data into pyef object
  datafile = datafolder+"/"+typ+nrun+".data.bin"
  if os.path.isfile(datafile) is False:
     print('File ',datafile,'does not exist... Wrong run type (Usage: >readData RUNID [RUNTYPE] [BOARDID])?')
     return

  print("Loading",datafile,"...")
  pyf = pyef.read_file(datafile)  #TBD: add error message if fails.
  print("Done.")
  return pyf


if __name__ == '__main__':
     if len(sys.argv)<2:
       print("Usage: >readData RUNID [RUNTYPE] [BOARDID]")
     else:

       #filename = "R{0}_trig_delays.npz".format(sys.argv[1])
       #fitDelays(filename)

       # First load data
       if len(sys.argv)>2:
         f = load_data(sys.argv[1],typ=sys.argv[2])
       else:
         f = load_data(sys.argv[1])
       if f == None:
         sys.exit()

       if len(sys.argv)==4:   # BOARDID parameter is present ==> display events
         print('Calling display_events() for unit',sys.argv[3])
         display_events(pyf = f,typ=sys.argv[2],tid=int(sys.argv[3]))
         sys.exit()

       # Perform coincidence search
       loadMaxCoarse(sys.argv[1])
       uid,distmat = build_distmat()
       trigtable = get_time(pyf=f)  # 2-lines matrix with [0,:]=UnitIDs and [1,:]=trigtimes
       #build_coincs(trigtable,uid,distmat)
