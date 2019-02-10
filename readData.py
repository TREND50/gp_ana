import sys
sys.path.append("../pyef/")
import pyef

import numpy as np
import pylab as pl

pl.ion()
c0 = 299792458
DISPLAY = 1
datafolder = "/home/martineau/GRAND/GRANDproto35/data/ulastai"

    
def twos_comp(val, bits):
    """compute the 2's compliment of int value val"""
    if (val & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
        val = val - (1 << bits)        # compute negative value
    return val   
    


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
  ntrigs = np.shape(trigtable)[1]
  tmax = np.max(d)/c0*1e9*10
  ants = trigtable[0,:]
  times = trigtable[1,:]
  i = 0
  while i<ntrigs: 
    trig_ref = times[i]
    ant_ref = ants[i]
    tsearch = times[i:-1]-times[i]
    tsearch = tsearch[np.where(tsearch<tmax)]  #  causal timewindow
    idsearch = ants[i:i+len(tsearch)]
    #print(tsearch,idsearch)
    tsearch = tsearch[np.where(idsearch!=ant_ref)]  # Remove triggers from target antenna
    if len(tsearch)>1: 
      # there are events in the causal timewindow
      print("Possible coinc!",tsearch,ants[i:i+len(tsearch)])
      i = i+1
    else:
      i = i+1
    #print(times[i])
    #print(trig_search)

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

  nsecs = []
  ttimes = []
  IDs = []
  for evt in f.event_list:
  # Loop on all events
    #print("\n\n!!!New event!!!")
    for ls in evt.local_station_list:
    # Loop on all units involved in event (should be one only at this stage)
      #ls.display()
      #print("LS ID=")
      #print(ls.header.ls_id)
      IDs.append(ls.header.ls_id) # 16 lowest bits of IP adress
    h = evt.header
    sec=h.event_sec
    if ls.header.ls_id==367:
      sec=sec+1  # Dirty fix R87!!! To be solved
    if ls.header.ls_id==387:
      sec=sec+16  # Dirty fix R91!!! To be solved
    nsec=h.event_nsec
    nsecs.append(nsec)
    ttime = sec*1e9+nsec
    ttimes.append(ttime)
    print("ID=",ls.header.ls_id,",Time=",sec,nsec)
    #evt.display()
    #h.display()
    
  nsecs = np.array(nsecs)
  ttimes = np.array(ttimes)
  ttimesns = ttimes-ttimes[0]
  ttimes = ttimes/1e9
  dur = max(ttimes)-min(ttimes)
  IDs = np.array(IDs)
  IDs = IDs-356  # --256 to go down to 8 lowest digits of IP adress & -100 to go down to antenna ID
  units = np.unique(IDs)
  res = np.vstack((IDs,ttimesns))
  res = np.sort(res)
  
  if DISPLAY:
    #pl.figure(1)
    #pl.subplot(211)
    #pl.plot(nsecs)
    #pl.subplot(212)
    #pl.hist(nsecs,100)
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

def load_data(nrun):
# Loads pyef object
  datafile = datafolder+"/R"+nrun+".data.bin"
  print("Loading",datafile,"...")
  pyf = pyef.read_file(datafile)  #TBD: add error message if fails.
  print("Done.")
  return pyf


if __name__ == '__main__':
     if len(sys.argv)!=2:
       print("Usage: >readData RUNID")
     else: 
       f = load_data(sys.argv[1])
       display_events(pyf = f,tid=9)
       uid,distmat = build_distmat()
       trigtable = get_time(pyf=f)  # 2-lines matrix with [0,:]=UnitIDs and [1,:]=trigtimes
       build_coincs(trigtable,uid,distmat)
     #pl.figure(1)
     #pl.plot(trigtable[1,:])
     #pl.figure(2)
     #pl.plot(trigtable[0,:])
     #pl.show()
