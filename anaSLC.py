# Script to analyze GP35 slow control data 
# OMH Aug. 29, 2018
# Updated Nov 16, on same structure as minBias analysis. Now analysing one board at a time only + reduce info to one txt file.

import os
from os.path import expanduser
import time
import datetime
#import mx.DateTime
import sys
import math
import numpy as np
import pylab as pl
import matplotlib.style
import matplotlib as mpl
import yaml
mpl.style.use('classic')


def loopSLCRuns(boardID,startrun,endrun):
  print("Calling loopRuns(). Will analyse minBias for board {0} between R{1} and R{2}.".format(boardID,startrun,endrun))
  time.sleep(1)
  
  for run in range(int(startrun),int(endrun)+1):
    loopSLCEvents(boardID,run)
    
def loopSLCEvents(boardID,RUNID):
   #datadir = "/home/pastsoft/data/"
   datadir = "/home/martineau/GRAND/GRANDproto35/data/ulastai/"
   filename = datadir+"S"+str(RUNID)+".yaml"   # To be modified
   if os.path.isfile(filename) is False:
     print('File ',filename,'does not exist. Aborting.')
     return

   resfile = 'SLC_b'+boardID+'.txt'  # Output file
   reso = open(resfile,'ab')
   a = np.loadtxt(resfile)
   a = np.loadtxt(resfile)
   try:
     tfmax =  a[-1,0]   
   except IndexError:  # When file is empty 
     tfmax = 0

   # Read data
   print('Scanning data file',filename)
   print("Loading data...")
   dataf=yaml.load_all(open(filename))
   print("Done.")
   
   IP = []
   VPower = []
   Th = []
   Temp = []
   maxCoarse = []
   TrigRate = []
   utcsec = []
   time_str = []
   for d in dataf:
     #determine whether it is a data:
     if d['msg_type']=='SLC':
       #print(d.keys())
       # Select data from this antenna only      
       uid=(d['source_ip'])[3]-100
       if uid != int(boardID):
	 #print 'This is board {0}, skiping it (analysing board {1} only)'.format(board,boardID)
         continue
       # Select unprocessed data only 
       thisUTC = d['received_timestamp'][0]
       if thisUTC<=tfmax: # Only looking at data more recent than already present in minBias_b[ID].txt
         print('Older data than in {0}, skiping it.'.format(resfile))
         continue	             
       if len(utcsec)>0 and thisUTC-utcsec[-1]<1:
         #print('Echoed data for unit {0}, skiping it.'.format(uid))
         continue	             

       utcsec.append(thisUTC)
       IP.append(d['source_ip'])
       power = [d['vpower1'], d['vpower2'],d['vpower3'], d['vpower4'],d['vpower5'], d['vpower6']]
       VPower.append(power)
       th = [d['th1m'], d['th1p'],d['th2m'], d['th2p'],d['th3m'], d['th3p']]
       Th.append(th)
       Temp.append(d['temp'])
       trate = [d['total_trig_rate'],d['ch1p_trig_rate'],d['ch1m_trig_rate'],d['ch2p_trig_rate'],d['ch2m_trig_rate'],d['ch3p_trig_rate'],d['ch3m_trig_rate']]
       TrigRate.append(trate)
       maxCoarse.append(d['max_coarse'])
       time_str.append(d['received_timestamp_str'])
       print(d['received_timestamp_str'])
             
	      
   # Write to file
   IP = np.asarray(IP,dtype=int)
   utcsec = np.asarray(utcsec,dtype=int)
   TrigRate = np.asarray(TrigRate,dtype=float)
   Th = np.asarray(Th,dtype=float)
   Temp = np.asarray(Temp,dtype=float)
   VPower = np.asarray(VPower,dtype=float)
   maxCoarse = np.asarray(maxCoarse,dtype=int)
   nev = len(utcsec)
   conc = np.concatenate((utcsec.reshape(nev,1),Temp.reshape(nev,1),VPower.reshape(nev,6),TrigRate.reshape(nev,7),maxCoarse.reshape(nev,1),),axis=1)   # Concatenate results
   conc = conc.reshape(np.size(utcsec),16) # 
   np.savetxt(reso,conc,fmt='%3.2f')  # Write to file

def displaySLC(boardID):
   home = expanduser("~")
   #resdir = home+"/GRAND/GRANDproto35/data/ulastai/"
   resdir = "./"
   resfile = resdir+"SLC_b"+str(boardID)+".txt"
   print("Calling displaySLC(). Will display SLC result file {0}".format(resfile))
   
   # Load data
   a = np.loadtxt(resfile)
   utc = a[:,0]
   temp = a[:,1]
   V = a[:,2:8]
   trig = a[:,8:16]
   mCoarse = a[:,15]
   sel = np.where(np.diff(utc)>0)  # remove redundant points
   utc = utc[sel[0]]
   temp = temp[sel[0]]
   V = V[sel[0],:]
   trig = trig[sel[0],:]
   mCoarse = mCoarse[sel[0]]
   
   sd,sm,sy=25,11,2018  # Start day,month,year
   ed,em,ey=18,11,2019  # End day,month,year
   startwindow=(datetime.datetime(sy,sm,sd)-datetime.datetime(1970,1,1)).total_seconds()
   endwindow=(datetime.datetime(ey,em,ed)-datetime.datetime(1970,1,1)).total_seconds()
   sel = np.where((utc<endwindow) & (utc> startwindow))
   utc = utc[sel]
   V = V[sel]
   temp = temp[sel]
   trig = trig[sel]
   mCoarse = mCoarse[sel]

   #time = (time-min(time[time>0]))/60
   Triglabel = ['Total','Ch1+','Ch2+','Ch3+','Ch1-','Ch2-','Ch3-']
   Voltlabel = ['Main','-3V','+4V','LNA1','LNA2','LNA3']
   print('SLC info for board',boardID,' in period:')
   print(datetime.datetime.fromtimestamp(min(utc)).strftime('%y/%m/%d - %H:%M:%S UTC'),' to ',datetime.datetime.fromtimestamp(max(utc)).strftime('%y/%m/%d - %H:%M:%S UTC'))
   print(np.shape(sel)[1], 'data points for board',boardID)
   datestart = datetime.datetime.fromtimestamp(min(utc)).strftime('%y/%m/%d %H:%M UTC')
   dateend = datetime.datetime.fromtimestamp(max(utc)).strftime('%y/%m/%d %H:%M UTC')
   print("Actual period displayed: {0}-{1}".format(datestart,dateend))

   # Time ticks
   nticks = 8
   ind = np.linspace(min(utc),max(utc),nticks)
   date = [datetime.datetime.fromtimestamp(ux).strftime('%H:%M') for ux in ind]
   #date = [datetime.datetime.fromtimestamp(ux).strftime('%m/%d') for ux in ind]
   
   # Temperature plot
   pl.figure(1)
   pl.plot(utc,temp)
   pl.xticks(np.linspace(min(utc),max(utc),nticks), date)
   pl.xlim(min(utc)-1,max(utc)+1)
   pl.grid(True)
   pl.legend(loc='best')
   pl.xlabel('Date [Month/Day]',size='large')
   #pl.xlabel('Time')
   pl.ylabel('Temperature ($^{\circ}$C)')
   pl.title('Board temperature')
   pl.savefig('temp.png')	

   # Voltage plot
   pl.figure(2)
   for i in range(np.shape(V)[1]):
     sub=321+i
     pl.subplot(sub)
     pl.plot(utc,V[:,i],lw=2,label=Voltlabel[i])
     pl.grid(True)
     if i>3:
       pl.xlabel('Date [Month/Day]',size='large')
       #pl.xlabel('Time')
     pl.xticks(np.linspace(min(utc),max(utc),nticks), date)
     pl.title(Voltlabel[i])
     pl.ylabel('Voltage (V)')
     pl.savefig('voltage.png')	

   pl.figure(3)  #Trig Rate
   # Plotting total trig rate only. 
   # pl.plot(time[sel],TrigRate[sel,0][0],lw=2,label=id) 
   # Plotting all individual trig rates. 
   labs = ['Total','X+','Y+','Z+','X-','Y-','Z-']
   for ch in range(7):
     pl.plot(utc,trig[:,ch],lw=2,label=labs[ch])

   pl.xticks(np.linspace(min(utc),max(utc),nticks), date)
   pl.grid(True)
   pl.ylabel('Total trig rate (Hz)')
   pl.xlabel('Time')
   #pl.xlabel('Date [Month/Day]',size='large')
   pl.legend(loc='best')
   pl.savefig('trig.png')	

   pl.draw()
   pl.show()

   return

def twos_comp(val, bits):
   """compute the 2's compliment of int value val"""
   if (val & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
     val = val - (1 << bits)	    # compute negative value
   return val


if __name__ == '__main__':
       loopSLCEvents(sys.argv[1],sys.argv[2])
       #loopSLCRuns(sys.argv[1],sys.argv[2],sys.argv[3])
       displaySLC(sys.argv[1])
