# -*- coding: utf-8 -*-
"""
Created on Mon May  8 09:57:44 2017

@author: hxu
"""
import numpy as np
from datetime import datetime,timedelta
import matplotlib.pyplot as plt
from math import sqrt
from netCDF4 import Dataset
def sh_rmtide(f,dt=1.,ends=0.):
    """
    removes solar and lunar tidal signal by sequentially averaging 
    over their periods: of 24h and 24h50.47m. This is equivalent to
    applying convolution with a trapezoidal shaped window (almost triagular).
    
    f - input tidal sequence
    dt - uniform time step of input series, in hours, dt=1. hr default
    ends - fill value for the ends: 0. or nan
    
    fm = sh_rmtide(f,dt=1.,ends=np.nan)
    """
    TS=24. # Solar hour angle period 24h
    TL=24.+50.47/60. # Lunar hour angle period
    N=len(f)
    fm=np.zeros(N)*ends

# remove solar period    
    T=TS/dt
    m=int(np.floor((T-1.)/2.))
    w=(T-(2*m+1))*0.5
# all weights =1 except at ends      
#    print w
    for i in range(m+1,N-m-1-1):
        #print len(f[i-m:i+m+1])
        fm[i]=(np.sum(f[i-m:i+m+1])+w*f[i-m-1]+w*f[i+m+1])/T
        
# remove lunar period
    f=fm*1.0  # deep copy!  
    T=TL/dt
    m=int(np.floor((T-1.)/2.))
    w=(T-(2*m+1))*0.5
# all weights =1 except at ends      
#    print w
    for i in range(m+1,N-m-1-1):
        #print len(f[i-m:i+m+1])
        fm[i]=(np.sum(f[i-m:i+m+1])+w*f[i-m-1]+w*f[i+m+1])/T
    
    return fm
URL='E0130.aanderaa.realtime.nc'
ds = Dataset(URL,'r').variables
f5=np.load('FVCOM_GOM3_GoMOOS_E01_201306.npy').tolist()

f6=np.load('FVCOM_GOM3_GoMOOS_E01_201307.npy').tolist()

f7=np.load('FVCOM_GOM3_GoMOOS_E01_201308.npy').tolist()

f8=np.load('FVCOM_GOM3_GoMOOS_E01_201309.npy').tolist()


f9=np.load('FVCOM_GOM3_GoMOOS_E01_201309.npy').tolist()

f10=np.load('FVCOM_GOM3_GoMOOS_E01_201310.npy').tolist()

f11=np.load('FVCOM_GOM3_GoMOOS_E01_201311.npy').tolist()


f12=np.load('FVCOM_GOM3_GoMOOS_E01_201312.npy').tolist()
f13=np.load('FVCOM_GOM3_GoMOOS_E01_201401.npy').tolist()

f14=np.load('FVCOM_GOM3_GoMOOS_E01_201402.npy').tolist()

f15=np.load('FVCOM_GOM3_GoMOOS_E01_201403.npy').tolist()
f16=np.load('FVCOM_GOM3_GoMOOS_E01_201404.npy').tolist()

f17=np.load('FVCOM_GOM3_GoMOOS_E01_201405.npy').tolist()

f18=np.load('FVCOM_GOM3_GoMOOS_E01_201406.npy').tolist()
'''
f19=np.load('FVCOM_GOM3_GoMOOS_E01_201306.npy').tolist()
'''
"""
plt.figure()
us=ds['current_u'][:,0,0,0]
ub=ds['current_u'][:,37,0,0]
ua=np.mean(ds['current_u'][:,:,0,0],1)

td=ds['time']
#units: days since 1858-11-17 00:00:00


plt.plot(t,us,'r-')
plt.plot(t,ub,'b-')
plt.plot(t,ua,'g-')
plt.show()

lat=ds['lat'];lon=ds['lon']
"""
t0=datetime(1858,11,17,0,0,0)
t=np.array([t0+timedelta(days=ds['time'][kt]) for kt in np.arange(len(ds['time']))])
index5=[]
index6=[]
index7=[]
index8=[]
index9=[]
index10=[]
u=[]
v=[]
time=[]
for a in np.arange(len(f5['t'])):
    if f5['t'][a]>=t[0] and f5['t'][a]<=t[-1]:
        u.append(f5['u0'][a])
        v.append(f5['v0'][a])
        time.append(f5['t'][a])
print len(u)

for a in np.arange(1,len(f6['t']),1):
    #if f6['t'][a]>=t[0] and f6['t'][a]<=t[-1]:
    u.append(f6['u0'][a])
    v.append(f6['v0'][a])
    time.append(f6['t'][a])
print len(u)

for a in np.arange(1,len(f7['t']),1):
    if f7['t'][a]>=t[0] and f7['t'][a]<=t[-1]:
        u.append(f7['u0'][a])
        v.append(f7['v0'][a])
        time.append(f7['t'][a])
print len(u)

for a in np.arange(1,len(f8['t']),1):
    if f8['t'][a]>=t[0] and f8['t'][a]<=t[-1]:
        u.append(f8['u0'][a])
        v.append(f8['v0'][a])
        time.append(f8['t'][a])


for a in np.arange(1,len(f9['t']),1):
    if f9['t'][a]>=t[0] and f9['t'][a]<=t[-1]:
        u.append(f9['u0'][a])
        v.append(f9['v0'][a])
        time.append(f9['t'][a])
print len(u)

for a in np.arange(1,len(f10['t']),1):
    if f10['t'][a]>=t[0] and f10['t'][a]<=t[-1]:
        u.append(f10['u0'][a])
        v.append(f10['v0'][a])
        time.append(f10['t'][a])

for a in np.arange(1,len(f11['t']),1):
    if f11['t'][a]>=t[0] and f11['t'][a]<=t[-1]:
        u.append(f11['u0'][a])
        v.append(f11['v0'][a])
        time.append(f11['t'][a])
print len(u)


for a in np.arange(1,len(f12['t']),1):
    if f12['t'][a]>=t[0] and f12['t'][a]<=t[-1]:
        u.append(f12['u0'][a])
        v.append(f12['v0'][a])
        time.append(f12['t'][a])
for a in np.arange(1,len(f13['t']),1):
    if f13['t'][a]>=t[0] and f13['t'][a]<=t[-1]:
        u.append(f13['u0'][a])
        v.append(f13['v0'][a])
        time.append(f13['t'][a])
print len(u)


for a in np.arange(1,len(f14['t']),1):
    if f14['t'][a]>=t[0] and f14['t'][a]<=t[-1]:
        u.append(f14['u0'][a])
        v.append(f14['v0'][a])
        time.append(f14['t'][a])

for a in np.arange(1,len(f15['t']),1):
    if f15['t'][a]>=t[0] and f15['t'][a]<=t[-1]:
        u.append(f15['u0'][a])
        v.append(f15['v0'][a])
        time.append(f15['t'][a])
for a in np.arange(1,len(f16['t']),1):
    if f16['t'][a]>=t[0] and f16['t'][a]<=t[-1]:
        u.append(f16['u0'][a])
        v.append(f16['v0'][a])
        time.append(f16['t'][a])
print len(u)


for a in np.arange(1,len(f17['t']),1):
    if f17['t'][a]>=t[0] and f17['t'][a]<=t[-1]:
        u.append(f17['u0'][a])
        v.append(f17['v0'][a])
        time.append(f17['t'][a])

for a in np.arange(1,len(f18['t']),1):
    if f18['t'][a]>=t[0] and f18['t'][a]<=t[-1]:
        u.append(f18['u0'][a])
        v.append(f18['v0'][a])
        time.append(f18['t'][a])
print len(u)

'''
for a in np.arange(1,len(f19['t']),1):
    if f19['t'][a]>=t[0] and f19['t'][a]<=t[-1]:
        u.append(f19['u0'][a])
        v.append(f19['v0'][a])
        time.append(f19['t'][a])
'''
currentu=[]
for a in np.arange(len(ds['current_u'])):
    currentu.append(ds['current_u'][a][0][0][0]/100)
currentv=[]
for a in np.arange(len(ds['current_u'])):
    currentv.append(ds['current_v'][a][0][0][0]/100)
fig,axes=plt.subplots(3,1,figsize=(7,10),sharex=True)
#plt.figure()

axes[0].set_title('u')
axes[0].plot(time,sh_rmtide(u),label='fvcom')
axes[0].plot(t,sh_rmtide(currentu),label='mooring E')
axes[0].set_ylabel('''m/s''')
axes[0].legend(loc='best')
#plt.show()
axes[1].set_title('v')
axes[1].plot(time,sh_rmtide(v),label='fvcom')
axes[1].plot(t,sh_rmtide(currentv),label='mooring E')
axes[1].set_ylabel('''m/s''')
axes[1].legend(loc='best')
speedf=[]
speedm=[]
for a in np.arange(len(sh_rmtide(v))):
    print 'a',a
    speedf.append(sqrt(sh_rmtide(v)[a]**2+sh_rmtide(u)[a]**2))
for a in np.arange(len(sh_rmtide(currentv))):
    print 'a',a
    speedm.append(sqrt(sh_rmtide(currentv)[a]**2+sh_rmtide(currentu)[a]**2))
axes[2].set_title('speed')
axes[2].set_ylabel('''m/s''')
axes[2].plot(time,speedf,label='fvcom')
axes[2].plot(t,speedm,label='mooring E')
axes[2].legend(loc='best')
plt.savefig('mooring E2013-6-24 TO 2014-6-2',dpi=400)

"""
axes[0].set_title('u')
axes[0].plot(time,u,label='fvcom')
axes[0].plot(t,currentu,label='mooring E')
axes[0].set_ylabel('''m/s''')
axes[0].legend(loc='best')
#plt.show()
axes[1].set_title('v')
axes[1].plot(time,v,label='fvcom')
axes[1].plot(t,currentv,label='mooring E')
axes[1].set_ylabel('''m/s''')
axes[1].legend(loc='best')
speedf=[]
speedm=[]
for a in np.arange(len(v)):
    print 'a',a
    speedf.append(sqrt(v[a]**2+u[a]**2))
for a in np.arange(len(currentv)):
    print 'a',a
    speedm.append(sqrt(currentv[a]**2+currentu[a]**2))
axes[2].set_title('speed')
axes[2].set_ylabel('''m/s''')
axes[2].plot(time,speedf,label='fvcom')
axes[2].plot(t,speedm,label='mooring E')
axes[2].legend(loc='best')
plt.savefig('mooring E2004-10-6 TO 2004-11-19',dpi=700)
"""