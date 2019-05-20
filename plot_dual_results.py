# takelot pfesi 36959
# File: plot_dual_resullts.py
#
# Desc: Plots the PSS of dual beam scans
#                                                                                                                                                                                                             
#=============================================================================#

#******************************************************************************
# Import packages/modules                                                     *
#******************************************************************************
#Built-in to python
#------------------
import numpy as np
import matplotlib.pylab as pl
import sys

#--------------
# User defined 
#--------------
from src.plot_functions import *
from src.exam_rms import *
from src.calc_pss import *

#--------------
# Code         
#--------------

#Create an argparse 
# Get source information from user

src = "3C123" #"3C123"
src2 = "3C123"
frq = '8400'

if frq=='5000':
 n='6'
elif frq=='8400':
 n='35'

# Setup file locations and output file
fileloc = 'datafiles/'+src+'/'+src2+'_'+n+'_plot_stats.csv'
dataloc = 'datafiles/'+src+'/'+src2+'_'+n+'_csvFiles/'+src2+'_'+n+'_filedata_db.csv'
outfile = src+'_'+n+'cm_DATA.csv'

# Get MSc data
#--------------
#data = np.genfromtxt(fileloc,dtype=str,delimiter=',',skip_header=2)
mdate=[]

mjdal=[]
mjdar=[]
mjdbl=[]
mjdbr=[]

mpssal=[]
mpssar=[]
mpssbl=[]
mpssbr=[]

ahpnl=[]
ahpsl=[]
aonl=[]
bhpnl=[]
bhpsl=[]
bonl=[]

ahpnr=[]
ahpsr=[]
aonr=[]

bhpnr=[]
bhpsr=[]
bonr=[]

mahpnl=[]
mahpsl=[]
maonl=[]

mbhpnl=[]
mbhpsl=[]
mbonl=[]

mahpnr=[]
mahpsr=[]
maonr=[]

mbhpnr=[]
mbhpsr=[]
mbonr=[]


#mdata =  'MSc_'+src+'_'+n+'_cm_PSS.csv'
mdata =  'MSc_HydraA_'+n+'_cm_PSS.csv'

with open(mdata,'r') as md:

    for _ in xrange(2):
        next(md)

    for line in md:
        m = line.split(',')[2]

        try: 
            alhpn= float((line.split(',')[25]).rstrip())
            mahpnl.append(float(m))
            ahpnl.append(float(alhpn))
        except Exception as e:
            pass

        try: 
            aslcp= float((line.split(',')[26]).rstrip())
            mahpsl.append(float(m))
            ahpsl.append(float(aslcp))
        except Exception as e:
            pass

        try: 
            aolcp= float((line.split(',')[27]).rstrip())
            maonl.append(float(m))
            aonl.append(float(aolcp))
        except Exception as e:
            pass
        
        try: 
            pal= float((line.split(',')[30]).rstrip())
            mjdal.append(float(m))
            mpssal.append(float(pal))
        except Exception as e:
            pass

        try: 
            bnlcp= float((line.split(',')[32]).rstrip())
            mbhpnl.append(float(m))
            bhpnl.append(float(bnlcp))
        except Exception as e:
            pass

        try:
            bslcp= float((line.split(',')[33]).rstrip())
            mbhpsl.append(float(m))
            bhpsl.append(float(bslcp))
        except:
            pass

        try:
            bolcp= float((line.split(',')[34]).rstrip())
            mbonl.append(float(m))
            bonl.append(float(bolcp))
        except:
            pass

        try:
            pbl = float((line.split(',')[37]).rstrip())
            mjdbl.append(float(m))
            mpssbl.append(float(pbl))
        except:
            pass

        try:
            anrcp= float((line.split(',')[39]).rstrip())
            ahpnr.append(float(m))
            mahpnr.append(float(anrcp))
        except:
            pass

        try:
            asrcp= float((line.split(',')[40]).rstrip())
            ahpsr.append(float(m))
            mahpsr.append(float(asrcp))
        except:
            pass

        try:
            aorcp= float((line.split(',')[41]).rstrip())
            aonr.append(float(m))
            maonr.append(float(aorcp))
        except:
            pass

        try:
            par = float((line.split(',')[44]).rstrip())
            mjdar.append(float(m))
            mpssar.append(float(par))
        except:
            pass

        try:
            bnrcp= float((line.split(',')[46]).rstrip())
            bhpnr.append(float(m))
            mbpnr.append(float(bnrcp))
        except:
            pass

        try:
            bsrcp= float((line.split(',')[47]).rstrip())
            bhpsr.append(float(m))
            mbhpsr.append(float(bsrcp))
        except:
            pass

        try:
            borcp= float((line.split(',')[48]).rstrip())
            bonr.append(float(m))
            mbonr.append(float(borcp))
        except:
            pass

        try:
            pbr = float((line.split(',')[51]).rstrip())
            mjdbr.append(float(m))
            mpssbr.append(float(pbr))
        except:
            pass

#sys.exit()
#mj=[]
#day0=[]
tl=[]
dtl=[]
tr=[]
dtr=[]
flux=[]
freq=[]
filen=[]

#Get data from the csv file
'''with open(dataloc,'r') as f:
 for _ in xrange(2):
  next(f)

 for line in f:

  m = line.split(',')[2]            # mjd
  d = line.split(',')[1]            # date
  fr=line.split(',')[4]             # frequency
  fl = float(line.split(',')[6])    # flux densitty, OTT

  day0.append(d)

  mj.append(float(m))'''


#print len(day0)
#sys.exit()
# Get data from plot file
hpnAlcp=[]
errhpnAlcp=[]
hpnBlcp=[]
errhpnBlcp =[]

hpnArcp=[]
errhpnArcp=[]
hpnBrcp=[]
errhpnBrcp =[]

hpsAlcp=[]
errhpsAlcp=[]
hpsBlcp=[]
errhpsBlcp =[]

hpsArcp=[]
errhpsArcp=[]
hpsBrcp=[]
errhpsBrcp =[]

onAlcp=[]
erronAlcp=[]
onBlcp=[]
erronBlcp =[]

onArcp=[]
erronArcp=[]
onBrcp=[]
erronBrcp =[]

mjd=[]
day=[]
ha=[]
pwv=[]

with open(fileloc,'r') as f:

 for _ in xrange(2):
  next(f)

 for line in f:

  fn     = line.split(',')[0]                 # filename
  filen.append(fn)

  d      = line.split(',')[1]                 # date
  day.append(d)

  m      = line.split(',')[2]                 # mjd
  mjd.append(float(m))

  fr     = line.split(',')[3]                 # frequency
  freq.append(fr)

  fl     = line.split(',')[4]                 # flux
  flux.append(float(fl))

  p      = line.split(',')[5]                 # pwv
  pwv.append(float(p))

  h      = line.split(',')[6]                 # ha
  ha.append(float(h))

  tal    = line.split(',')[7]                 # tsys lcp  
  tl.append(float(tal))

  dtal   = (line.split(',')[8]).strip()       # dtsys lcp
  dtl.append(float(dtal))

  tar    = line.split(',')[9]                 #tsys rcp
  tr.append(float(tar))

  dtar   = (line.split(',')[10]).strip()      #dtsys rcp
  dtr.append(float(dtar))

  hpnal  = line.split(',')[11]                # hpn A beam lcp
  ehpnal = (line.split(',')[12]).strip()      # err hpn A beam lcp
  hpnbl  = line.split(',')[13]                # hpn B beam lcp
  ehpnbl = (line.split(',')[14]).strip()      # err hpn B beam lcp

  hpnAlcp.append(float(hpnal))
  errhpnAlcp.append(float(ehpnal))
  hpnBlcp.append(float(hpnbl))
  errhpnBlcp.append(float(ehpnbl))

  hpnar  = line.split(',')[20]                # hpn A beam rcp
  ehpnar = (line.split(',')[21]).strip()      # err hpn A beam rcp
  hpnbr  = line.split(',')[22]                # hpn B beam rcp
  ehpnbr = (line.split(',')[23]).strip()      # err hpn B beam rcp

  hpnArcp.append(float(hpnar))
  errhpnArcp.append(float(ehpnar))
  hpnBrcp.append(float(hpnbr))
  errhpnBrcp.append(float(ehpnbr))

  hpsal  = line.split(',')[29]
  ehpsal = (line.split(',')[30]).strip()
  hpsbl  = line.split(',')[31]
  ehpsbl = (line.split(',')[32]).strip() 

  hpsAlcp.append(float(hpsal))
  errhpsAlcp.append(float(ehpsal))
  hpsBlcp.append(float(hpsbl))
  errhpsBlcp.append(float(ehpsbl))

  hpsar  = line.split(',')[38]
  ehpsar = (line.split(',')[39]).strip()
  hpsbr  = line.split(',')[40]
  ehpsbr = (line.split(',')[41]).strip() 

  hpsArcp.append(float(hpsar))
  errhpsArcp.append(float(ehpsar))
  hpsBrcp.append(float(hpsbr))
  errhpsBrcp.append(float(ehpsbr))

  onal  = line.split(',')[47]
  eonal = (line.split(',')[48]).strip()
  onbl  = line.split(',')[49]
  eonbl = (line.split(',')[50]).strip() 

  onAlcp.append(float(onal))
  erronAlcp.append(float(eonal))
  onBlcp.append(float(onbl))
  erronBlcp.append(float(eonbl))

  onar  = line.split(',')[56]
  eonar = (line.split(',')[57]).strip()
  onbr  = line.split(',')[58]
  eonbr = (line.split(',')[59]).strip() 

  onArcp.append(float(onar))
  erronArcp.append(float(eonar))
  onBrcp.append(float(onbr))
  erronBrcp.append(float(eonbr))

#print len(mj), len(mjd)
#print mj[0], mjd[0]
#sys.exit()
#plot beams
'''pl.plot(mbonl,bonl,'b.')
pl.plot(mjd,onBlcp,'g.')
pl.savefig('Beams.eps')
pl.show()
pl.close()
#print len(mj),len(mjd)
sys.exit()'''

# CALCULATE PSS
pal=[]
epal=[]
pbl=[]
epbl=[]
par=[]
epar=[]
pbr=[]
epbr=[]

tcoral=[]
tcorar=[]
tcorbl=[]
tcorbr=[]
tcorale=[]
tcorare=[]
tcorble=[]
tcorbre=[]
pcalcp=[]
pcarcp=[]
pcblcp=[]
pcbrcp=[]

for i in range(len(hpnAlcp)):
 #print hpsAlcp[i],hpsBlcp[i],hpsArcp[i],hpsBrcp[i]
 # PSS LCP
 pssAlcp, psserrAlcp, pcAlcp, TcorrAlcp, TcorerrAlcp = cals_pss_db(hpnAlcp[i], errhpnAlcp[i], hpsAlcp[i], errhpsAlcp[i], onAlcp[i], erronAlcp[i],fl)
 pssBlcp, psserrBlcp, pcBlcp, TcorrBlcp, TcorerrBlcp = cals_pss_db(hpnBlcp[i], errhpnBlcp[i], hpsBlcp[i], errhpsBlcp[i], onBlcp[i], erronBlcp[i],fl)
# PSS RCP
 pssArcp, psserrArcp, pcArcp, TcorrArcp, TcorerrArcp = cals_pss_db(hpnArcp[i], errhpnArcp[i], hpsArcp[i], errhpsArcp[i], onArcp[i], erronArcp[i],fl)
 pssBrcp, psserrBrcp, pcBrcp, TcorrBrcp, TcorerrBrcp = cals_pss_db(hpnBrcp[i], errhpnBrcp[i], hpsBrcp[i], errhpsBrcp[i], onBrcp[i], erronBrcp[i],fl)

 #print pssAlcp
 #sys.exit()
 pal.append(pssAlcp)
 epal.append(psserrAlcp)
 pbl.append(pssBlcp)
 epbl.append(psserrBlcp )
 par.append(pssArcp)
 epar.append(psserrArcp)
 pbr.append(pssBrcp)
 epbr.append(psserrBrcp)
 tcoral.append(TcorrAlcp)
 tcorar.append(TcorrArcp)
 tcorbl.append(TcorrBlcp)
 tcorbr.append(TcorrBrcp)
 tcorale.append(TcorerrAlcp)
 tcorare.append(TcorerrArcp)
 tcorble.append(TcorerrBlcp)
 tcorbre.append(TcorerrBrcp)
 pcalcp.append(pcAlcp)
 pcarcp.append(pcArcp)
 pcblcp.append(pcBlcp)
 pcbrcp.append(pcBrcp)

#print len(mj),len(mjd)
#print len(mjd), len(pal)
#sys.exit()

#---------------------------------------------------------
# Calculate mean values, remove outliers
#-----------------------------------------------------------
mpsal,psal,tal,sdminal,sdmaxal,meanal = exam(mjd,pal,len(pal))
mpsar,psar,tar,sdminar,sdmaxar,meanar = exam(mjd,par,len(par))
mpsbl,psbl,tbl,sdminbl,sdmaxbl,meanbl = exam(mjd,pbl,len(pbl))
mpsbr,psbr,tbr,sdminbr,sdmaxbr,meanbr = exam(mjd,pbr,len(pbr))

#sys.exit()

dt1=[]
dt2=[]
dt3=[]
dt4=[]

date1=[]
date2=[]
date3=[]
date4=[]

# Get dates for all pss
for i in range(len(mjd)):
 for j in  range(len(mpsal)):
  if(mjd[i]==mpsal[j]):
   dt1.append(epal[i])
   date1.append(day[i])

for i in range(len(mjd)):
 for j in  range(len(mpsar)):
  if(mjd[i]==mpsar[j]):
    dt3.append(epar[i])
    date3.append(day[i])

for i in range(len(mjd)):
 for j in  range(len(mpsbl)):
  if(mjd[i]==mpsbl[j]):
   dt2.append(epbl[i])
   date2.append(day[i])

for i in range(len(mjd)):
 for j in  range(len(mpsbr)):
  if(mjd[i]==mpsbr[j]):
    dt4.append(epbr[i])
    date4.append(day[i])

#print len(dt4), len(mpsbr)
#sys.exit()

#Create arrays to store values for standard deviation cut-off points
mn1= np.ones(len(mpsal))*meanal
sm1= np.ones(len(mpsal))*sdminal
sx1= np.ones(len(mpsal))*sdmaxal
mn3= np.ones(len(mpsar))*meanar
sm3= np.ones(len(mpsar))*sdminar
sx3= np.ones(len(mpsar))*sdmaxar
mn2= np.ones(len(mpsbl))*meanbl
sm2= np.ones(len(mpsbl))*sdminbl
sx2= np.ones(len(mpsbl))*sdmaxbl
mn4= np.ones(len(mpsbr))*meanbr
sm4= np.ones(len(mpsbr))*sdminbr
sx4= np.ones(len(mpsbr))*sdmaxbr

# PLOT PSS RESULTS
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12
pl.rc('font', size=SMALL_SIZE)          # controls default text sizes
pl.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
pl.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
pl.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
pl.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
pl.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
pl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
#pl.rc('title', size=BIGGER_SIZE)          # controls default text sizes
figure()

pl.subplot(411)
pl.title('PSS plots Automated vs Manual data ', fontsize=BIGGER_SIZE)
pl.ylabel('A_BEAM_LCP  [Jy/K]')
pl.plot(mpsal,mn1,'k',label="Mean PSS = %.1f" %(meanal))
pl.plot(mpsal,sm1,'gray')
pl.plot(mpsal,sx1,'gray')
pl.legend(loc="best")
pl.plot(mpsal,psal,'b.')
#pl.errorbar(mpsal,psal,yerr=dt1,fmt="b.")
#pl.plot(mjdal,mpssal,'k.')

pl.subplot(412)
pl.xlabel('Time [MJD]')
pl.ylabel('B_BEAM_LCP  [Jy/K]')
pl.plot(mpsbl,mn2,'k',label="Mean PSS = %.1f" %(meanbl))
pl.plot(mpsbl,sm2,'gray')
pl.plot(mpsbl,sx2,'gray')
pl.legend(loc="best")
pl.plot(mpsbl,psbl,'g.')
#pl.errorbar(mpsbl,psbl,yerr=dt2,fmt="g.")
#pl.plot(mjdbl,mpssbl,'k.')

pl.subplot(413)
pl.ylabel('A_BEAM_RCP  [Jy/K]')
pl.plot(mpsar,mn3,'k',label="Mean PSS = %.1f" %(meanar))
pl.plot(mpsar,sm3,'gray')
pl.plot(mpsar,sx3,'gray')
pl.legend(loc="best")
pl.plot(mpsar,psar,'r.')
#pl.errorbar(mpsar,psar,yerr=dt3,fmt="r.")
#pl.plot(mjdar,mpssar,'k.')

pl.subplot(414)
pl.xlabel('Time [MJD]')
pl.ylabel('B_BEAM_RCP  [Jy/K]')
pl.plot(mpsbr,mn4,'k',label="Mean PSS = %.1f" %(meanbr))
pl.plot(mpsbr,sm4,'gray')
pl.plot(mpsbr,sx4,'gray')
pl.legend(loc="best")
#pl.errorbar(mpsbr,psbr,yerr=dt4,fmt="m.")
pl.plot(mpsbr,psbr,'m.')
#pl.plot(mjdbr,mpssbr,'k.')
pl.savefig(src+'_'+n+'cm_PSS_PLOTS.eps')
#pl.show()
pl.close()

# WRITE FINAL OUTPUT TO FILE
with open(outfile,'w') as f:
    f.write('%s,%s,%s,%s,%s,%s,%s,%s,\
    %s,%s,%s,%s,\
    %s,%s,%s,%s,%s,%s,\
    %s,%s,%s,%s,%s,\
    %s,%s,%s,%s,%s,%s,\
    %s,%s,%s,%s,%s,\
    %s,%s,%s,%s,%s,%s,\
    %s,%s,%s,%s,%s,\
    %s,%s,%s,%s,%s,%s,\
    %s,%s,%s,%s,%s\n' 
    %('Filename'  , 'Date'  , 'MJD'  , 'SRC'   , 'FREQ'     , 'FLUX'       , 'PWV', 'HA',
    'TsysLCP'        , 'dTsysLCP'          , 'TsysRCP'        , 'dTsysRCP'              , 
    'A_BEAM_HPN_LCP' , 'errA_BEAM_HPN_LCP' , 'A_BEAM_HPS_LCP' , 'errA_BEAM_HPS_LCP'     , 'A_BEAM_ON_LCP'               , 'errA_BEAM_ON_LCP' , 
    'PC'             , 'corrON_LCP'        , 'err_corrON_LCP' , 'A_BEAM_PSS_LCP [Jy/K]' , 'error A_BEAM_PSS_LCP [Jy/K]' ,

    'B_BEAM_HPN_LCP' , 'errB_BEAM_HPN_LCP' , 'B_BEAM_HPS_LCP' , 'errB_BEAM_HPS_LCP'     , 'B_BEAM_ON_LCP'               , 'errB_BEAM_ON_LCP' , 
    'PC'             , 'corrON_LCP'        , 'err_corrON_LCP' , 'B_BEAM_PSS_LCP [Jy/K]' , 'error B_BEAM_PSS_LCP [Jy/K]' ,

    'A_BEAM_HPN_RCP' , 'errA_BEAM_HPN_RCP' , 'A_BEAM_HPS_RCP' , 'errA_BEAM_HPS_RCP'     , 'A_BEAM_ON_RCP'               , 'errA_BEAM_ON_RCP' , 
    'PC'             , 'corrON_RCP'        , 'err_corrON_RCP' , 'A_BEAM_PSS_RCP [Jy/K]' , 'error A_BEAM_PSS_RCP [Jy/K]' ,

    'B_BEAM_HPN_RCP' , 'errB_BEAM_HPN_RCP' , 'B_BEAM_HPS_RCP' , 'errB_BEAM_HPS_RCP'     , 'B_BEAM_ON_RCP'               , 'errB_BEAM_ON_RCP' , 
    'PC'             , 'corrON_RCP'        , 'err_corrON_RCP' , 'B_BEAM_PSS_RCP [Jy/K]' , 'error B_BEAM_PSS_RCP [Jy/K]' ,

))

    day=sorted(day)
    #print date1#date1=sorted(date1)
    #sys.exit()
    for j in range (len(day)):

        a=0
        b=0
        c=0
        z=0

        for i in range(len(date1)):
            if(day[j]==date1[i]):
                a = date1[i]
                break

        for p in range(len(date2)):
            if day[j]==date2[p] :
                b = date2[p] 
                break

        for t in range(len(date3)):
            if day[j]==date3[t]:
                c = date3[t]
                break

        for g in range(len(date4)):
            if day[j]==date4[g]:
                z = date4[g]
                break

        #print '\n-->> ', 'j: ',j,', day[j]: ',day[j],', len(j): ', len(day),  ', a: ',a,  ', b: ',b ,  ', c: ',c,  ', z: ',z #freq[j],src#,'\n'#0,mjd[j]
        #sys.exit()
        if(day[j] ==a and day[j]!=b and day[j]==c and day[j]==z):

             print 'A0CZ ', a,b ,c,z

             f.write('%s,%s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s\n'
                    %(filen[j] , day[j]        , mjd[j]     , src           , freq[j], 
                    flux[0]    , pwv[j], ha[j], tl[j]         , dtl[j]     , tr[j]         , dtr[j], 
                    hpnAlcp[j] , errhpnAlcp[j] , hpsAlcp[j] , errhpsAlcp[j] , onAlcp[j], erronAlcp[j],  
                    pcalcp[j]  , tcoral[j]     , tcorale[j] , psal[i]       , dt1[i]   , 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    hpnArcp[j] , errhpnArcp[j] , hpsArcp[j] , errhpsArcp[j] , onArcp[j], erronArcp[j],  
                    pcarcp[j]  ,tcorar[j]      , tcorare[j] , psar[t]       , dt3[t]   , 
                    hpnBrcp[j] , errhpnBrcp[j] , hpsBrcp[j] , errhpsBrcp[j] , onBrcp[j], erronBrcp[j], 
                    pcbrcp[j]  ,tcorbr[j]      , tcorbre[j] , psbr[g]       ,dt4[g]
                    ))

        elif(day[j] == a and day[j] == b and day[j] != c and day[j] == z):

           print' AB0Z - ',a,b,c,z

           f.write('%s,%s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s\n'
                    %(filen[j] , day[j]        , mjd[j]     , src           , freq[j], 
                    flux[0]    , pwv[j], ha[j], tl[j]         , dtl[j]     , tr[j]         , dtr[j], 
                    hpnAlcp[j] , errhpnAlcp[j] , hpsAlcp[j] , errhpsAlcp[j] , onAlcp[j], erronAlcp[j],  
                    pcalcp[j]  , tcoral[j]     , tcorale[j] , psal[i]       , dt1[i]   , 
                    hpnBlcp[j] , errhpnBlcp[j] , hpsBlcp[j] , errhpsBlcp[j] , onBrcp[j], erronBlcp[j],  
                    pcblcp[j]  ,tcorbl[j]      , tcorble[j] , psbl[p]       , dt2[p]   , 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    hpnBrcp[j] , errhpnBrcp[j] , hpsBrcp[j] , errhpsBrcp[j] , onBrcp[j], erronBrcp[j], 
                    pcbrcp[j]  ,tcorbr[j]      , tcorbre[j] , psbr[g]       ,dt4[g]
                    ))

        elif(day[j] == a and day[j] == b and day[j] == c and day[j] != z):

            '''
            print ' ABC0 - ', a,b,c,z
            print 'file: ',filen[j]
            print 'day: ',day[j]
            print 'mjd: ', mjd[j]
            print 'hpnalcp: ',  hpnAlcp[j] 
            print 'err: ', errhpnAlcp[j] 
            print 'hpsalcp: ', hpsAlcp[j] 
            print 'err: ', errhpsAlcp[j] 
            print 'onalcp: ', onAlcp[j] 
            print 'err: ',erronAlcp[j]
            print 'pcAlcp: ', pcalcp[j]  
            print 'tcoral: ', tcoral[j]     
            print 'err: ', tcorale[j] 
            print 'psal: ', psal[i]       
            print 'dt1: ', dt1[i]   
            print 'hpnblcp: ', hpnBlcp[j] 
            print 'err: ', errhpnBlcp[j] 
            print 'hpnblcp: ', hpsBlcp[j] 
            print 'err: ', errhpsBlcp[j] 
            print 'onbrcp: ', onBrcp[j]
            print 'err: ', erronBlcp[j]
            print 'pcblcp: ', pcblcp[j]  
            print 'tcorbl: ',tcorbl[j]      
            print 'tcorble: ', tcorble[j] 
            print 'psbl: ', psbl[p]
            print 'dt3: ' , dt2[p]   
            print 'hpnarcp: ',  hpnArcp[j] 
            print 'err: ', errhpnArcp[j]  
            print 'hpsarcp: ',hpsArcp[j] 
            print 'err: ', errhpsArcp[j] 
            print 'onalcp: ', onArcp[j]
            print 'err: ', erronArcp[j]
            print 'pcarcp: ', pcarcp[j]  
            print 'tcorar: ',tcorar[j]      
            print 'tcorare: ', tcorare[j] 
            print 'psar: ', psar[t]      
            print 'dt2: ', dt3[t]   
            print 'hpnbrcp: ', hpnBrcp[j] 
            print 'err: ', errhpnBrcp[j] 
            print 'hpsbrcp: ', hpsBrcp[j] 
            print 'err: ', errhpsBrcp[j] 
            print 'onbrcp: ', onBrcp[j]
            print 'err: ', erronBrcp[j]'''
            f.write('%s,%s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s\n'
                    %(filen[j] , day[j]        , mjd[j]     , src   , freq[j], 
                    flux[0]    , pwv[j], ha[j], tl[j]         , dtl[j]     , tr[j]         , dtr[j], 
                    hpnAlcp[j] , errhpnAlcp[j] , hpsAlcp[j] , errhpsAlcp[j] , onAlcp[j], erronAlcp[j],  
                    pcalcp[j]  , tcoral[j]     , tcorale[j] , psal[i]       , dt1[i]   , 
                    hpnBlcp[j] , errhpnBlcp[j] , hpsBlcp[j] , errhpsBlcp[j] , onBrcp[j], erronBlcp[j],  
                    pcblcp[j]  ,tcorbl[j]      , tcorble[j] , psbl[p]     , dt2[p]   , 
                    hpnArcp[j] , errhpnArcp[j] , hpsArcp[j] , errhpsArcp[j] , onArcp[j], erronArcp[j],  
                    pcarcp[j]  ,tcorar[j]      , tcorare[j] , psar[t]      , dt3[t]   , 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', ))
 
        elif(day[j] == a and day[j] == b and day[j] != c and day[j] != z):

            print ' AB00 - ',a,b,c,z      

            f.write('%s,%s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s\n'
                    %(filen[j] , day[j]        , mjd[j]     , src           , freq[j], 
                    flux[0]    , pwv[j], ha[j], tl[j]         , dtl[j]     , tr[j]         , dtr[j], 
                    hpnAlcp[j] , errhpnAlcp[j] , hpsAlcp[j] , errhpsAlcp[j] , onAlcp[j], erronAlcp[j],  
                    pcalcp[j]  , tcoral[j]     , tcorale[j] , psal[i]       , dt1[i]   , 
                    hpnBlcp[j] , errhpnBlcp[j] , hpsBlcp[j] , errhpsBlcp[j] , onBrcp[j], erronBlcp[j],  
                    pcblcp[j]  ,tcorbl[j]      , tcorble[j] , psbl[p]       , dt2[p]   , 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    ))

        elif(day[j] == a and day[j] == b and day[j] == c and day[j] == z):

            print ' ABCZ - ',a,b,c,z

            f.write('%s,%s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s\n'
                    %(filen[j] , day[j]        , mjd[j]     , src           , freq[j], 
                    flux[0]    ,pwv[j], ha[j],  tl[j]         , dtl[j]     , tr[j]         , dtr[j], 
                    hpnAlcp[j] , errhpnAlcp[j] , hpsAlcp[j] , errhpsAlcp[j] , onAlcp[j], erronAlcp[j],  
                    pcalcp[j]  , tcoral[j]     , tcorale[j] , psal[i]       , dt1[i]   , 
                    hpnBlcp[j] , errhpnBlcp[j] , hpsBlcp[j] , errhpsBlcp[j] , onBrcp[j], erronBlcp[j],  
                    pcblcp[j]  ,tcorbl[j]      , tcorble[j] , psbl[p]       , dt2[p]   , 
                    hpnArcp[j] , errhpnArcp[j] , hpsArcp[j] , errhpsArcp[j] , onArcp[j], erronArcp[j],  
                    pcarcp[j]  ,tcorar[j]      , tcorare[j] , psar[t]       , dt3[t]   , 
                    hpnBrcp[j] , errhpnBrcp[j] , hpsBrcp[j] , errhpsBrcp[j] , onBrcp[j], erronBrcp[j], 
                    pcbrcp[j]  ,tcorbr[j]      , tcorbre[j] , psbr[g]       ,dt4[g]
                    ))

        elif(day[j] == a and day[j] != b and day[j] != c and day[j] == z):

            print ' A00Z - ',a,b,c,z

            f.write('%s,%s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s\n'
                    %(filen[j] , day[j]        , mjd[j]     , src           , freq[j], 
                    flux[0]    , pwv[j], ha[j], tl[j]         , dtl[j]     , tr[j]         , dtr[j], 
                    hpnAlcp[j] , errhpnAlcp[j] , hpsAlcp[j] , errhpsAlcp[j] , onAlcp[j], erronAlcp[j],  
                    pcalcp[j]  , tcoral[j]     , tcorale[j] , psal[i]       , dt1[i]   , 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    hpnBrcp[j] , errhpnBrcp[j] , hpsBrcp[j] , errhpsBrcp[j] , onBrcp[j], erronBrcp[j], 
                    pcbrcp[j]  ,tcorbr[j]      , tcorbre[j] , psbr[g]       ,dt4[g]
                    ))

        elif(day[j] == a and day[j] != b and day[j] == c and day[j] != z):

            print ' A0C0 - ',a,b,c,z

            f.write('%s,%s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s\n'
                    %(filen[j] , day[j]        , mjd[j]     , src           , freq[j], 
                    flux[0]    , pwv[j], ha[j], tl[j]         , dtl[j]     , tr[j]         , dtr[j], 
                    hpnAlcp[j] , errhpnAlcp[j] , hpsAlcp[j] , errhpsAlcp[j] , onAlcp[j], erronAlcp[j],  
                    pcalcp[j]  , tcoral[j]     , tcorale[j] , psal[i]       , dt1[i]   , 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    hpnArcp[j] , errhpnArcp[j] , hpsArcp[j] , errhpsArcp[j] , onArcp[j], erronArcp[j],  
                    pcarcp[j]  ,tcorar[j]      , tcorare[j] , psar[t]       , dt3[t]   , 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-',   
                    ))

        elif(day[j] == a and day[j] != b and day[j] != c and day[j] != z):

            print ' A000 - ',a,b,c,z

            f.write('%s,%s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s\n'
                    %(filen[j] , day[j]        , mjd[j]     , src           , freq[j], 
                    flux[0]    , pwv[j], ha[j], tl[j]         , dtl[j]     , tr[j]         , dtr[j], 
                    hpnAlcp[j] , errhpnAlcp[j] , hpsAlcp[j] , errhpsAlcp[j] , onAlcp[j], erronAlcp[j],  
                    pcalcp[j]  , tcoral[j]     , tcorale[j] , psal[i]       , dt1[i]   , 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-',    
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-',                  ))

        elif(day[j] != a and day[j] == b and day[j] == c and day[j] == z):

            print ' 0BCZ - ',a,b,c,z

            f.write('%s,%s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s\n'
                    %(filen[j] , day[j]        , mjd[j]     , src           , freq[j], 
                    flux[0]    , pwv[j], ha[j], tl[j]         , dtl[j]     , tr[j]         , dtr[j], 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    hpnBlcp[j] , errhpnBlcp[j] , hpsBlcp[j] , errhpsBlcp[j] , onBrcp[j], erronBlcp[j],  
                    pcblcp[j]  ,tcorbl[j]      , tcorble[j] , psbl[p]       , dt2[p]   , 
                    hpnArcp[j] , errhpnArcp[j] , hpsArcp[j] , errhpsArcp[j] , onArcp[j], erronArcp[j],  
                    pcarcp[j]  ,tcorar[j]      , tcorare[j] , psar[t]       , dt3[t]   , 
                    hpnBrcp[j] , errhpnBrcp[j] , hpsBrcp[j] , errhpsBrcp[j] , onBrcp[j], erronBrcp[j], 
                    pcbrcp[j]  ,tcorbr[j]      , tcorbre[j] , psbr[g]       ,dt4[g]
                    ))

        elif(day[j] != a and day[j] == b and day[j] != c and day[j] == z):

            print ' 0B0Z - ',a,b,c,z

            f.write('%s,%s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s\n'
                    %(filen[j] , day[j]        , mjd[j]     , src           , freq[j], 
                    flux[0]    , pwv[j], ha[j], tl[j]         , dtl[j]     , tr[j]         , dtr[j], 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    hpnBlcp[j] , errhpnBlcp[j] , hpsBlcp[j] , errhpsBlcp[j] , onBrcp[j], erronBlcp[j],  
                    pcblcp[j]  ,tcorbl[j]      , tcorble[j] , psbl[p]       , dt2[p]   , 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    hpnBrcp[j] , errhpnBrcp[j] , hpsBrcp[j] , errhpsBrcp[j] , onBrcp[j], erronBrcp[j], 
                    pcbrcp[j]  ,tcorbr[j]      , tcorbre[j] , psbr[g]       ,dt4[g]
                    ))

        elif(day[j] != a and day[j] != b and day[j] == c and day[j] == z):

            print ' 00CZ - ',a,b,c,z

            f.write('%s,%s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s\n'
                    %(filen[j] , day[j]        , mjd[j]     , src           , freq[j], 
                    flux[0]    , pwv[j], ha[j], tl[j]         , dtl[j]     , tr[j]         , dtr[j], 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    hpnBlcp[j] , errhpnBlcp[j] , hpsBlcp[j] , errhpsBlcp[j] , onBrcp[j], erronBlcp[j],  
                    '-'  , '-'  , '-','-' , '-'    , 
                    hpnArcp[j] , errhpnArcp[j] , hpsArcp[j] , errhpsArcp[j] , onArcp[j], erronArcp[j],  
                    pcarcp[j]  ,tcorar[j]      , tcorare[j] , psar[t]       , dt3[t]   , 
                    hpnBrcp[j] , errhpnBrcp[j] , hpsBrcp[j] , errhpsBrcp[j] , onBrcp[j], erronBrcp[j], 
                    pcbrcp[j]  ,tcorbr[j]      , tcorbre[j] , psbr[g]       ,dt4[g]
                    ))

        elif(day[j] != a and day[j] != b and day[j] != c and day[j] == z):

            print ' 000Z - ',a,b,c,z

            f.write('%s,%s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s\n'
                    %(filen[j] , day[j]        , mjd[j]     , src           , freq[j], 
                    flux[0]    ,pwv[j], ha[j],  tl[j]         , dtl[j]     , tr[j]         , dtr[j], 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    hpnBrcp[j] , errhpnBrcp[j] , hpsBrcp[j] , errhpsBrcp[j] , onBrcp[j], erronBrcp[j], 
                    pcbrcp[j]  ,tcorbr[j]      , tcorbre[j] , psbr[g]       ,dt4[g]
                    ))

        elif(day[j] != a and day[j] == b and day[j] == c and day[j] != z):

            print ' 0BC0 ',a,b,c,z

            f.write('%s,%s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s\n'
                    %(filen[j] , day[j]        , mjd[j]     , src           , freq[j], 
                    flux[0]    , pwv[j], ha[j], tl[j]         , dtl[j]     , tr[j]         , dtr[j], 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    hpnBlcp[j] , errhpnBlcp[j] , hpsBlcp[j] , errhpsBlcp[j] , onBrcp[j], erronBlcp[j],  
                    pcblcp[j]  ,tcorbl[j]      , tcorble[j] , psbl[p]       , dt2[p]   , 
                    hpnArcp[j] , errhpnArcp[j] , hpsArcp[j] , errhpsArcp[j] , onArcp[j], erronArcp[j],  
                    pcarcp[j]  ,tcorar[j]      , tcorare[j] , psar[t]       , dt3[t]   , 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    ))

        elif(day[j] != a and day[j] != b and day[j] == c and day[j] != z):

            print ' 00C0 - ',a,b,c,z

            f.write('%s,%s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s\n'
                    %(filen[j] , day[j]        , mjd[j]     , src           , freq[j], 
                    flux[0]    ,pwv[j], ha[j], tl[j]         , dtl[j]     , tr[j]         , dtr[j], 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    hpnArcp[j] , errhpnArcp[j] , hpsArcp[j] , errhpsArcp[j] , onArcp[j], erronArcp[j],  
                    pcarcp[j]  ,tcorar[j]      , tcorare[j] , psar[t]       , dt3[t]   , 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    ))

        elif(day[j] != a and day[j] == b and day[j] != c and day[j] != z):

            print ' 0B00 - ',a,b,c,z

            f.write('%s,%s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,\
                    %s,%s,%s,%s,%s,%s\n'
                    %(filen[j] , day[j]        , mjd[j]     , src           , freq[j], 
                    flux[0]    , pwv[j], ha[j], tl[j]         , dtl[j]     , tr[j]         , dtr[j], 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    hpnBlcp[j] , errhpnBlcp[j] , hpsBlcp[j] , errhpsBlcp[j] , onBrcp[j], erronBlcp[j],  
                    pcblcp[j]  ,tcorbl[j]      , tcorble[j] , psbl[p]       , dt2[p]   , 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    '-'  , '-'  , '-','-' , '-',  '-',
                    '-'  , '-'  , '-','-' , '-', 
                    ))

        elif (day[j] != b and day[j] !=c and day[j] !=z and day[j] !=a):
            print '# no match'
            pass

        else:
            print "I don't understand you daddy! "

    f.close()

