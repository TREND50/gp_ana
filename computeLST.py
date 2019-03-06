import time
import ephem
import datetime
import sys

# Compute LST time for given location

print("Call is: >python computeLST.py year month day [h min s]")

# Gansu 
longi = "94.10"
lat = "39.35"
#Ulastai
#longi = "86.71"
#lat = "42.95"
# Nancay
longi = "2.11"
lat = "47.2"
y = int(sys.argv[1])
m = int(sys.argv[2])
d = int(sys.argv[3])
h, mn, s = 8, 0 , 0  # UTC<=>GMT time
if len(sys.argv)>4:
  h = int(sys.argv[4])
  if len(sys.argv)>5:
    mn = int(sys.argv[5])
    if len(sys.argv)>6:
      s = int(sys.argv[6])
  
site = ephem.Observer();
site.date = datetime.datetime(y,m,d,h,mn,s)   # UTC<=>GMT time
site.long = ephem.degrees(longi)
site.lat = ephem.degrees(lat)
site.elevation = 1500;

print("Site: ({0}N,{1}E)".format(longi,lat))
print("UTC time: {0}/{1}/{2} {3}:{4}:{5}".format(y,m,d,h,mn,s))
print("Local Sideral time:",site.sidereal_time())
