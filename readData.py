import os
import sys
sys.path.append("../")
import pyef

import numpy as np
import pylab as pl
import yaml

pl.ion()
c0 = 299792458
DISPLAY = 0
datafolder = "/home/martineau/GRAND/GRANDproto35/data/ulastai"
#datafolder = "/mnt/disk"
#IDsin = []
utcSLC = []    
maxCoarse = []
def twos_comp(val, bits):
    """compute the 2's compliment of int value val"""
    if (val & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
        val = val - (1 << bits)        # compute negative value
    return val   


def loadMaxCoarse(runid): 
  global IDsin
  allIDs = []
  allUTC = []
  allMaxCoarse = []
  slcfile = datafolder+"/S"+runid+".yaml"
  if os.path.isfile(slcfile) is False:
     print('File ',slcfile,'does not exist. Aborting.')
     return
  
  # Now read SLC file ad dump infos in arrays
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
  #try:
  #  if len(IDsin) == 0:  # SLC data was not yet loaded
  #    print("Ooops... No MaxCoarse info yet! Fetching it from SLC data.")
  #    loadMaxCoarse(sys.argv[1])     
  #  else:
  #    print("MaxCoarse info available for following units",IDsin)
     
  #except:
  #  print("Ooops... No MaxCoarse info yet! Fetching it from SLC data.")
  #  loadMaxCoarse(sys.argv[1])      
  
  
  # Now retrieve proper maxCoarse info
  i = np.nonzero(IDsin == uid)[0]  
  i = i[0]
  indt = np.argmin(np.abs(utcSLC[i]-utcsec))
  #print(utcsec,utcSLC[i][indt],maxCoarse[i][indt])
  return maxCoarse[i][indt]

  #print(utcSLC,maxCoarse)
  
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
  #p = [x, y ,z]  # Why cannot numpy build matrixes with this syntax?????????? Makes me crazy
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
  ntrigs = np.shape(trigtable)[0]
  tmax = np.max(d)/c0*1e9*1.5
  uid = trigtable[:,0]  # Vector of unit UDs
  times = trigtable[:,1]  # Vector of ttrig times (units = ns; ref = 1st trigger)
  i = 0
  #
  delays = []
  uids = []
  while i<ntrigs:   # Loop on all triggers in table
    trig_ref = times[i]
    id_ref = uid[i]
    tsearch = times[i:-1]-trig_ref
    tsearch = tsearch[np.argwhere(tsearch<tmax)].T[0]  #  causal timewindow. Transpose needed to get a line vector and avoid []
    idsearch = uid[i:i+len(tsearch)]
    #print("Units in causal timewindow:",idsearch,tsearch,times[i:i+len(tsearch)])
    others = np.argwhere(idsearch!=id_ref).T[0]  # Remove triggers from target antenna
    #print("From different units:",idsearch[others],tsearch[others],others)
    #print(tsearch,len(tsearch),np.size(tsearch),np.shape(tsearch))
    if len(tsearch)>1 and sum(tsearch)>0: 
      # there are events in the causal timewindow
      print("*** Reference unit:",id_ref,trig_ref,", now looking for coincs within",tmax,"ns")
      print("Possible coinc with units",idsearch[others],tsearch[others],others)
      uids.append(idsearch[others][0])
      delays.append(tsearch[others][0])
      i = i+others[-1]+1
    else:
      i = i+1
      
    #print(times[i])
    #print(trig_search)
  uids = np.array(uids)
  delays = np.array(delays)
  pl.figure(1)
  pl.hist(delays[uids==5],100)
  pl.hist(delays[uids==9],100)
  pl.hist(delays[uids==18],100)
  pl.show()
                    
def get_time(nrun=None,pyf=None):
# Retrieves time info from datafile
  
  if pyf == None:
    print("No pyef object, loading it from run number.")
    if nrun == None:
      print("get_time error! Pass run number or pyef object as argument")
      return
    pyf = load_data(nrun)

  
  nevts = len(pyf.event_list)
  #print(nevts,"events in file",datafile)
  # TBD: access file header

  secs = []
  nsecs = []
  ttimes = []
  IDs = []
  for evt in f.event_list:
  # Loop on all events
    #print("\n\n!!!New event!!!")
    for ls in evt.local_station_list:
    # Loop on all units involved in event (should be one only at this stage)
      #print("Local unit info")
      #ls.display()
      #print("LS ID=")
      #print(ls.header.ls_id)
      uid = int(ls.header.ls_id-356) # 16 lowest bits of IP adress --256 to go down to 8 lowest digits of IP adress & -100 to go down to antenna ID
      IDs.append(uid) 
      nsec = ls.header.gps_nanoseconds  # What is that one for???
      
    h = evt.header
    sec = h.event_sec
    if 0:
      if uid == 9:
        sec = sec-3  # Dirty fix R167!!! To be solved
      if uid == 18:
        sec = sec-5  # Dirty fix R167!!! To be solved
    nsec = h.event_nsec
    # Now correct time from maxCoarse value
    maxcoarse = getMaxCoarse(uid,sec)
    cor=125e6/(getMaxCoarse(uid,sec)+1)
    if abs(cor-1)>0.01: # Abnormal correction value
      cor = 1
    nsec = nsec*cor
    nsecs.append(nsec)
    ttime = sec*1e9+nsec
    ttimes.append(ttime)
    secs.append(sec)
    print("ID=",uid,",Time=",sec,nsec,h.event_nsec,nsec-h.event_nsec)
    #print("Event info")
    #evt.display()
    #print("Header info")
    #h.display()
    
  secs = np.array(secs)
  nsecs = np.array(nsecs)
  ttimes = np.array(ttimes)
  ttimesns = ttimes-min(ttimes)
  ttimes = ttimesns/1e9
  dur = max(ttimes)-min(ttimes)
  IDs = np.array(IDs)
  units = np.unique(IDs)
  ind = np.argsort(ttimesns)
  ttimesns_ordered = ttimesns[ind]
  IDs_ordered = IDs[ind]
  res = np.vstack((IDs_ordered,ttimesns_ordered))
  res = res.T
  #for a in range(1,np.shape(res)[0]):
  #  print(res[a,0],res[a,1],res[a,1]-res[a-1,1])
  for uid in units:
      tdif = np.diff(ttimes[IDs==uid])
      aid = np.argwhere(tdif<0)
      for i in aid:
        #print(ind)
        #i = ind[1]
        #print(i,j,aid[j],xx)
        print("*** Error for unit",uid,":\nevent",i,": SSS =",secs[IDs==uid][i],"\nevent",i+1,": SSS =",secs[IDs==uid][i+1])
	
  if DISPLAY:
    # Build delta_ns info
    for uid in units:
      deltat = []
      thisID = np.argwhere(IDs==uid)
      otherID = np.argwhere(IDs!=uid)
      for i in thisID:
        i = i[0]
        tg = min(abs(otherID-i))[0]  # Closest index of other units
        tg = tg+i
        print("****",i,tg)
        print(otherID[max(0,tg-10):tg+10])
        otherNS = nsecs[otherID[max(0,tg-10):tg+10]]
        print(otherNS)
        if len(otherNS)>0:
          dt = min(abs(otherNS-nsecs[i]))
          deltat.append(dt)
      
      pl.figure(12)
      pl.hist(deltat)
      pl.show()
      	
    pl.figure(1)
    pl.subplot(211)
    for uid in units:
      pl.plot(nsecs[IDs==uid],label=uid)
    pl.xlabel('Trigger nb')
    pl.ylabel('Nanosec counter value')
    pl.legend(loc='best')
    pl.subplot(212)
    pl.hist(nsecs,100)
    pl.xlabel('Nanosec counter value')
    pl.xlim([0,1e9])
    pl.figure(2)
    for uid in units:
      pl.plot(ttimes[IDs==uid],label=uid)
    pl.xlabel('Trigger nb')
    pl.ylabel('Trigger time (s)')
    pl.legend(loc='best')
    pl.title('Triggers')
    
    pl.figure(3)
    pl.hist(IDs,100)
    pl.xlabel("Unit ID")
    
    pl.show()

  print(nevts,"events in run",nrun)
  print("Run duration:",dur,"seconds.")
  print("Units present in run:")
  for uid in units:
    print("Unit",uid,":",np.shape(np.where(IDs==uid))[1],"events.")
  
  return res


def display_events(nrun=None,pyf=None,tid=None):
# Display events
  lab = ['X','Y','Z','Cal']
  if pyf == None:
    print("No pyef object, loading it from run number.")
    if nrun == None:
      print("get_time error! Pass run number or pyef object as argument")
      return  
    pyf = load_data(nrun)
  
  nevts = len(pyf.event_list)
  #print(nevts,"events in file",datafile)
  # TBD: access file header
  for evt in f.event_list:
  # Loop on all events
    #print("\n\n!!!New event!!!")
    for ls in evt.local_station_list:
    # Loop on all units involved in event (should be one only at this stage)
      uid = ls.header.ls_id - 356 # Remove 255 for digits 8-16 in IP adress & 100 for unit ID
      if tid==None or uid ==tid:  
      # Display 
        raw = ls.adc_buffer
        hraw = [hex(int(a)) for a in raw]  # Transfer back to hexadecimal
        draw = [twos_comp(int(a,16), 12) for a in hraw] #2s complements
        draw = np.array(draw)*1./2048  # in Volts
        nsamples = int(ls.header.trace_length/4)  # draw corresponds to 4 channels
        #offset = int(nsamples/2.0)  # Offset position at center of waveform
        #print nsamples,"samples per channel --> offset = ",offset
        thisEvent = np.reshape(draw,(4,nsamples));
        thisEvent = pow(10,(thisEvent+np.min(thisEvent)))
        tmus = np.array(range(nsamples))*20e-3  # Time axis in mus
        evtnb = ls.header.event_nr
        pl.figure(1)
        for i in range(3):
          pl.plot(tmus[3:],thisEvent[i][3:],label=lab[i])
        pl.plot([tmus[int(nsamples/2)+15], tmus[int(nsamples/2)+15]],[np.min(thisEvent[:][:]),np.max(thisEvent[:][:])])
        pl.title('Evt {0} Antenna {1}'.format(evtnb,uid))
        pl.xlim(tmus[3],max(tmus))
        pl.xlabel('Time ($\mu$s)')
        pl.ylabel('Voltage (V)')
        pl.legend(loc="best")
        pl.show()
        input()
        pl.close(1)
	
	
def load_data(nrun):
# Loads pyef object
  datafile = datafolder+"/R"+nrun+".data.bin"
  if os.path.isfile(datafile) is False:
     print('File ',datafile,'does not exist. Aborting.')
     return

  print("Loading",datafile,"...")
  pyf = pyef.read_file(datafile)  #TBD: add error message if fails.
  print("Done.")
  return pyf


if __name__ == '__main__':
     if len(sys.argv)!=2:
       print("Usage: >readData RUNID")
     else: 
       loadMaxCoarse(sys.argv[1])
       f = load_data(sys.argv[1])
       if f == None:
         sys.exit()
       #display_events(pyf = f,tid=9)
       uid,distmat = build_distmat()
       #print(uid,distmat)
       #input()
       trigtable = get_time(pyf=f)  # 2-lines matrix with [0,:]=UnitIDs and [1,:]=trigtimes
       build_coincs(trigtable,uid,distmat)
     #pl.figure(1)
     #pl.plot(trigtable[1,:])
     #pl.figure(2)
     #pl.plot(trigtable[0,:])
     #pl.show()
