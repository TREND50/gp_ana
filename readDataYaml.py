# Script to analyze GP35 data saved under yaml format (debug) 
# OMH March 5, 2019
# Taken from anaSLC.py

import os
from os.path import expanduser
import time
import datetime
import sys
import math
import numpy as np
import pylab as pl
import matplotlib.style
import matplotlib as mpl
import yaml
mpl.style.use('classic')
from dateutil.parser import parse
from readData import twos_comp

DISPLAY = 0
    
def loopEvents(RUNID):
   lab = ['X','Y','Z','Cal']
  
   #datadir = "/home/pastsoft/data/"
   datadir = "/home/martineau/GRAND/GRANDproto35/data/ulastai/"
   filename = datadir+"R"+str(RUNID)+".data.yaml"   
   if os.path.isfile(filename) is False:
     print('File ',filename,'does not exist. Aborting.')
     return

   # Read data
   print('Scanning data file',filename)
   print("Loading data...")
   dataf=yaml.load_all(open(filename))
   print("Done.")
   
   IP = []
   timetag = [];
   timens = [];
   times = []
   utcsec = []
   i = 0
   for d in dataf:
     #determine whether it is a data:
     if d['msg_type']=='DATA':
       #print(d.keys())
       # Select data from this antenna only      
       uid=(d['source_ip'])[3]-100
       evtid = d['event_count']
       thisUTC = d['received_timestamp'][0]
       utcsec.append(thisUTC)
       IP.append(d['source_ip'])
       sss = d['sss']
       times.append(sss) # trigger second info ("second since start")
       ts2 = d['ts2']
       ts1trig = d['ts1trigger']
       ts1pps = d['ts1pps']
       tns = (ts2*4+ts1pps-ts1trig)*2  # trigger ns info
       timens.append(tns)
       ttag = sss*1e9+tns  # triger time
       timetag.append(ttag) 
       
       sss_cor = d['sss_corr']
       timestr = d['received_timestamp_str']
       dateobj = parse(timestr)
       utc = dateobj.timestamp()
       
       print("Triger",i,"on unit",uid,"at time",timestr,sss,sss_cor)

       if DISPLAY: # Display events
         raw = d['data']
         hraw = [hex(int(a)) for a in raw]  # Transfer back to hexadecimal
         draw = [twos_comp(int(a,16), 12) for a in hraw] #2s complements
         draw = np.array(draw)*1./2048  # in Volts
         nsamples = 180  # harcoded
         thisEvent = np.reshape(draw,(4,nsamples));
         thisEvent = pow(10,(thisEvent+np.min(thisEvent)))
         tmus = np.array(range(nsamples))*20e-3  # Time axis in mus
         pl.figure(1)
         for j in range(3):
           pl.plot(tmus[3:],thisEvent[j][3:],label=lab[j])
         pl.plot([tmus[int(nsamples/2)+15], tmus[int(nsamples/2)+15]],[np.min(thisEvent[:][:]),np.max(thisEvent[:][:])])
         pl.title('Evt {0} Antenna {1}'.format(evtid,uid))
         pl.xlim(tmus[3],max(tmus))
         pl.xlabel('Time ($\mu$s)')
         pl.ylabel('Voltage (V)')
         pl.legend(loc="best")
         pl.show()
         input()
         pl.close(1)
		    
       i = i+1
       if i>=1000: # ùacx number of events to be read
         break
	     
	      
   # Write to file
   IP = np.asarray(IP,dtype=int)
   utcsec = np.asarray(utcsec,dtype=int)
   timetag = np.asarray(timetag)
   times = np.asarray(times)
   timens = np.asarray(timens)
   
   # Now display final time info
   pl.figure(1)
   pl.subplot(2,1,1)
   pl.plot(times)
   pl.xlabel('Trigger count')
   pl.ylabel('Second info')
   pl.subplot(2,1,2)
   pl.hist(timens)
   pl.xlabel('Nanosecond info')
   pl.show()
   
   

if __name__ == '__main__':
       loopEvents(sys.argv[1])
