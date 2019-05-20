#plot results - single beam
import numpy as np
import matplotlib.pylab as pl
from src.plot_functions import *
from src.exam_rms import *
#plot pss

src= "HydraA"
fileloc = 'datafiles/'+src+'/'+src+'_13NB_plot_stats.csv'
dataloc = 'datafiles/'+src+'/'+src+'_13NB_csvFiles/'+src+'_13NB_filedata_sb.csv'
outdata = src+'_13NB_DATA.csv'

#data = np.genfromtxt(fileloc,dtype=str,delimiter=',',skip_header=2)
#mdata =  np.genfromtxt(mscdatadtype=str,delimiter=',',skip_header=2)

mj=[]
tl=[]
dtl=[]
tr=[]
dtr=[]

with open(dataloc,'r') as f:
 for _ in xrange(10):
  next(f)
 for line in f:
  m = line.split(',')[2]
  tal =line.split(',')[19]
  dtal =(line.split(',')[20]).strip()
  tar =line.split(',')[21]
  dtar =(line.split(',')[22]).strip() 

  mj.append(float(tal))
  tl.append(float(dtal))
  dtl.append(float(tar))
  tr.append(float(dtar))
  dtr.append(float(m))

l=[]
r=[]
pssl=[]
pssr=[]
mjd=[]

with open(fileloc,'r') as f:
 for _ in xrange(2):
  next(f)
 for line in f:
  m = line.split(',')[2]
  tal =line.split(',')[5]
  dtal =(line.split(',')[6]).strip()
  tar =line.split(',')[14]
  dtar =(line.split(',')[15]).strip() 

  l.append(float(tal))
  pssl.append(float(dtal))
  r.append(float(tar))
  pssr.append(float(dtar))
  mjd.append(float(m))

mpl=np.mean(l)
mpr=np.mean(r)

#print l
t,t1,t2,sdmin,sdmax,mean = exam(mjd,l,len(l))
ts,t1s,t2s,sdmins,sdmaxs,means = exam(mjd,r,len(r))

#print t2,t2s
#sys.exit()


dt1=[]
dt2=[]

for i in range(len(mjd)):
 for j in  range(len(t)):
  if(mjd[i]==t[j]):
   dt1.append(pssl[i])

for i in range(len(mjd)):
 for j in  range(len(ts)):
  if(mjd[i]==ts[j]):
    dt2.append(pssr[i])

mn1= np.ones(len(t))*mean
sm1= np.ones(len(t))*sdmin
sx1= np.ones(len(t))*sdmax
mn2= np.ones(len(ts))*means
sm2= np.ones(len(ts))*sdmins
sx2= np.ones(len(ts))*sdmaxs

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12
pl.rc('font', size=SMALL_SIZE)          # controls default text sizes
pl.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
pl.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
pl.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
pl.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
pl.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
pl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#print len(ml),len(mpl)
figure()

pl.subplot(211)
pl.title('PSS plots')
pl.ylabel('PSS_LCP  [Jy/K]')
pl.plot(t,mn1,'k',label="Mean = %.1f" %(mean))
pl.plot(t,sm1,'r')
pl.plot(t,sx1,'r')
pl.legend(loc="best")
#pl.plot(mjd,l,'r.')
pl.errorbar(t,t1,yerr=dt1,fmt="b.")
pl.plot(t,t1,'b.')

'''pl.plot(t,mn1,'k',label="Mean = %.1f" %(mean))
pl.plot(t,sm1,'r')
pl.plot(t,sx1,'r')
pl.plot(ml,mpl,'k.',label="Msc data mean = %.1f" %(mld))
pl.legend(loc="best")
#pl.plot(x,pssrcp,'c.',label="PSS_RCP")
#pl.errorbar(x,pssrcp,yerr=pssrcperr,fmt="c.")
#pl.plot(xfr,yfr,'k.')'''

pl.subplot(212)
pl.xlabel('Time [MJD]')
pl.ylabel('PSS_RCP  [Jy/K]')
pl.plot(ts,mn2,'k',label="Mean = %.1f" %(means))
pl.plot(ts,sm2,'r')
pl.plot(ts,sx2,'r')
#pl.plot(ml,mpl,'k.',label="Msc data mean = %.1f" %(mld))
pl.legend(loc="best")
#pl.plot(mjd,r,'r.')
pl.errorbar(ts,t1s,yerr=dt2,fmt="g.")
pl.plot(ts,t1s,'g.')
pl.savefig(src+'_13NB_PSS_PLOTS.jpg')
pl.show()
'''pl.plot(ts,mn2,'k',label="Mean = %.1f" %(means))
pl.plot(ts,sm2,'r')
pl.plot(ts,sx2,'r')
pl.plot(mr,mpr,'k.',label="Msc data mean = %.1f" %(mrd))
pl.legend(loc="best")
pl.legend(loc='best')
pl.savefig('PSSmins')
#pl.show()
pl.close()


#Write out PSS value
with open("PSS_VALUES.csv",'w') as f:
 f.write('%s,%s,%s\n' %('Freq','PSSLCP','PSSRCP'))
 f.write('%s,%f,%f\n' %('2280',mean,means))
 
 
 #Other statistics
figure()
pl.title("Other statictics")

pl.subplot(421)
pl.title("Other statictics LCP")
pl.plot(mjd,splb,'b.',label="before fit")
pl.plot(mjd,spla,'g.',label="after fit")
pl.ylabel("Spline of data")
pl.legend(loc="best")

pl.subplot(423)
pl.plot(mjd,spld,'k.',label="difference")
pl.ylabel("Spline of data")
pl.legend(loc="best")

pl.subplot(425)
pl.plot(mjd,s2n,'k.')
pl.ylabel("S/N")

pl.subplot(427)
pl.plot(mjd,flag,'k.')
pl.ylabel("Flag")

pl.subplot(422)
pl.title("Other statictics RCP")
pl.plot(mjd,splbr,'b.',label="before fit")
pl.plot(mjd,splar,'g.',label="after fit")
pl.ylabel("Spline of data")
pl.legend(loc="best")

pl.subplot(424)
pl.plot(mjd,spldr,'k.',label="difference")
pl.ylabel("Spline of data")
pl.legend(loc="best")

pl.subplot(426)
pl.plot(mjd,s2nr,'k.')
pl.ylabel("S/N")

pl.subplot(428) 
pl.plot(mjd,flagr,'k.')
pl.ylabel("Flag")
pl.savefig("Otherstatistics.jpg")
#pl.show()
pl.close()
'''
