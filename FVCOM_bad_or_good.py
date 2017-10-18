# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 13:13:15 2017

@author: bling
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import  interpolate
from datetime import datetime, timedelta
URL1='current_04hind_hourly.nc'
t1=datetime(2001,7,1,0,0,0,0)
t2=datetime(2002,11,30,0,0,0,0)
################################
t3=datetime(2002,12,1,0,0,0,0)
t4=datetime(2003,3,31,0,0,0,0)
####################################
t5=datetime(2003,4,1,0,0,0,0)
t6=datetime(2009,3,31,0,0,0,0)
####################################
t7=datetime(2009,4,1,0,0,0,0)
t8=datetime(2010,7,31,0,0,0,0)

####################################
t9=datetime(2010,8,1,0,0,0,0)
t10=datetime(2011,6,30,0,0,0,0)
####################################
t11=datetime(2011,7,1,0,0,0,0)
t12=datetime(2014,5,31,0,0,0,0)
####################################
t13=datetime(2014,6,1,0,0,0,0)
t14=datetime(2014,10,31,0,0,0,0)
####################################
t15=datetime(2014,11,1,0,0,0,0)
t16=datetime(2015,5,30,0,0,0,0)
####################################
t17=datetime(2015,6,1,0,0,0,0)
t18=datetime(2015,11,30,0,0,0,0)
#####################################################################3
#########################################################################
######################################################################
t110=datetime(2001,7,1,0,0,0,0)
t21=datetime(2008,12,31,0,0,0,0)
################################
t31=datetime(2009,1,1,0,0,0,0)
t41=datetime(2009,3,31,0,0,0,0)
####################################
t51=datetime(2009,4,1,0,0,0,0)
t61=datetime(2009,5,31,0,0,0,0)
####################################
t71=datetime(2009,6,1,0,0,0,0)
t81=datetime(2010,5,31,0,0,0,0)

####################################
t91=datetime(2010,6,1,0,0,0,0)
t101=datetime(2011,1,31,0,0,0,0)
####################################
t111=datetime(2011,2,1,0,0,0,0)
t121=datetime(2012,10,31,0,0,0,0)
####################################
t131=datetime(2012,11,1,0,0,0,0)
t141=datetime(2012,12,31,0,0,0,0)
####################################
t151=datetime(2013,1,1,0,0,0,0)
t161=datetime(2013,2,28,0,0,0,0)
####################################
t171=datetime(2013,3,1,0,0,0,0)
t181=datetime(2013,6,30,0,0,0,0)

t191=datetime(2013,7,1,0,0,0,0)
t201=datetime(2013,11,30,0,0,0,0)

t211=datetime(2013,12,1,0,0,0,0)
t221=datetime(2015,1,31,0,0,0,0)

t231=datetime(2015,2,1,0,0,0,0)
t241=datetime(2015,12,31,0,0,0,0)

fig,axes=plt.subplots(1,1,figsize=(10,5))
plt.plot([t1,t2],[1,1],'r')

plt.plot([t3,t4],[1,1],'b',linewidth=4)

plt.plot([t5,t6],[1,1],'r')

plt.plot([t7,t8],[1,1],'b',linewidth=4)

plt.plot([t9,t10],[1,1],'r')

plt.plot([t11,t12],[1,1],'b',linewidth=4)

plt.plot([t13,t14],[1,1],'r')
plt.plot([t15,t16],[1,1],'b',linewidth=4)
plt.plot([t17,t18],[1,1],'r')


plt.plot([t110,t21],[0.5,0.5],'r',label='good')
plt.plot([t31,t41],[0.5,0.5],'b',linewidth=4,label='bad')
plt.plot([t51,t61],[0.5,0.5],'r')
plt.plot([t71,t81],[0.5,0.5],'b',linewidth=4)
plt.plot([t91,t101],[0.5,0.5],'r')
plt.plot([t111,t121],[0.5,0.5],'b',linewidth=4)
plt.plot([t131,t141],[0.5,0.5],'r')
plt.plot([t151,t161],[0.5,0.5],'b',linewidth=4)
plt.plot([t171,t181],[0.5,0.5],'r')
plt.plot([t191,t201],[0.5,0.5],'b',linewidth=4)
plt.plot([t211,t221],[0.5,0.5],'r')
plt.plot([t231,t241],[0.5,0.5],'b',linewidth=4)
plt.text(datetime(2009,1,1,0,0,0,0),1.2,'Mooring E')
plt.text(datetime(2009,1,1,0,0,0,0),0.7,'Mooring I')
plt.legend()
plt.ylim([0,1.5])
plt.title('Mooring VS FVCOM (2m)')
plt.savefig('mooring e',dpi=200)
