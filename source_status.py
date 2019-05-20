#P V van Zyl
#
# source_status.py
#****************************************************************************************************
#	IMPORT PACKAGES                                                                                                                                        *
#****************************************************************************************************

#Built-in to python
#------------------------
import numpy as np
import os
import sys

#Declare variable/s
CSV_TAG = 'csvFiles'

#Create output file/s
STATUS_FILE = 'source_status.csv'

#Set file location/s
STATUS_LOC = 'source_status/'

#create if not exists
with open(STATUS_LOC + STATUS_FILE , 'w') as f:
    f.write( '%s, %s, %s, %s, %s\n'%('Source', 'Freq', 'Start date', 'End date', 'No. Files'))

#Fetch source files
files = os.walk('files').next()

folders=files[1]
folders=sorted(folders)

fopen = open(STATUS_LOC+STATUS_FILE ,'a')

for i in range(len(folders)):

    subfolders = os.walk('files/'+folders[i]).next()
    for j in range(len(subfolders[1])):

        subfolders1 = sorted(subfolders[1])
      
        if CSV_TAG in subfolders1[j]:
            subfiles=  subfolders1[j]
            #subfiles=sorted(subfiles)
            p = os.walk('files/'+folders[i]+'/'+subfiles).next()[2]
            #print len(p),subfiles
            for k in range(len(p)):
                #print p[k]
                if 'filedata' in p[k]:
                    #print 'source: ', folders[i], ', folder: ', subfiles.strip('_csvFiles'),', file: ',p[k]

                    if folders[i]=='HydraA':

                        loc= 'files/'+folders[i]+'/'+subfiles+'/'+p[k]

                        if '13' in subfiles:
                            freq=2280
                            f=np.genfromtxt(loc,delimiter=',',dtype=str,skip_header=10)

                        elif ('6' in subfiles) :
                            freq=5000
                            f=np.genfromtxt(loc,delimiter=',',dtype=str,skip_header=2)

                        elif ('35' in subfiles):
                            freq=8400
                            f=np.genfromtxt(loc,delimiter=',',dtype=str,skip_header=2)
    
                        elif '2.5' in subfiles:
                            freq=12178
                            f=np.genfromtxt(loc,delimiter=',',dtype=str,skip_header=6)

                        loc= 'files/'+folders[i]+'/'+subfiles+'/'+p[k]

                        print '\nSource: ', folders[i], ', Folder: ', subfiles.strip('_csvFiles'), 'Freq: ',freq,', Start date: ',(f[0][1]), ', End date: ', f[-1][1], ', No Files: ',len(f),'\n'
                        fopen.write( '%s,%s,%s,%s,%s\n'%( folders[i],freq,(f[0][1]),f[-1][1],len(f)))

                    elif folders[i]=='Jupiter':
                        loc= 'files/'+folders[i]+'/'+subfiles+'/'+p[k]
                        f=np.genfromtxt(loc,delimiter=',',dtype=str,skip_header=23)
                        
                        freq=22000

                        print '\nSource: ', folders[i], ', Folder: ', subfiles.strip('_csvFiles'), 'Freq: ',freq,', Start date: ',(f[0][1]), ', End date: ', f[-1][1], ', No Files: ',len(f),'\n'
                        fopen.write( '%s,%s,%s,%s,%s\n'%( folders[i],freq,(f[0][1]),f[-1][1],len(f)))
                    else:
                        loc= 'files/'+folders[i]+'/'+subfiles+'/'+p[k]
                        if ('2280' in subfiles): 
                            freq=2280
                            f=np.genfromtxt(loc,delimiter=',',dtype=str,skip_header=6)

                        elif ('22000' in subfiles):
                            freq=22000
                            f=np.genfromtxt(loc,delimiter=',',dtype=str,skip_header=6)

                        elif ('5000' in subfiles): 
                            freq=5000
                            f=np.genfromtxt(loc,delimiter=',',dtype=str,skip_header=2)

                        elif ('8400' in subfiles):
                            freq=8400
                            f=np.genfromtxt(loc,delimiter=',',dtype=str,skip_header=2)

                        elif ('12178' in subfiles):
                            freq=12178
                            f=np.genfromtxt(loc,delimiter=',',dtype=str,skip_header=9)

                        if (len(f)==13):
                            fopen.write( '%s,%s,%s,%s,%s\n'%( folders[i],freq,(f[1]),f[1],1))
                            print '\nSource: ', folders[i], ', Folder: ', subfiles.strip('_csvFiles'), 'Freq: ',freq,', Start date: ',(f[1]), ', End date: ', f[1], ', No Files: ',1,'\n'                            
                        else:
                            fopen.write( '%s,%s,%s,%s,%s\n'%( folders[i],freq,(f[0][1]),f[-1][1],len(f)))
                            print '\nSource: ', folders[i], ', Folder: ', subfiles.strip('_csvFiles'), 'Freq: ',freq,', Start date: ',(f[0][1]), ', End date: ', f[-1][1], ', No Files: ',len(f),'\n'


