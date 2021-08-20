import csv
import numpy as np
import cv2
import sys
import math

fw=open("ratio.dat","w+")

opm_list=open("/localdata/zzhang624/data/ompx/PDBs_OPM/ModPDBsv2/files_name",'r')
content = opm_list.read()
pdbs = content.split("\n")
pdbs.pop()
print(pdbs)
opm_list.close()

for x in pdbs:
    fr=open("xy_dat/"+x+"_xy.dat",'r')
    xx=list(csv.reader(fr,delimiter=' ',quoting=csv.QUOTE_NONNUMERIC))
    yy=[]
    for xxx in xx:
      yyy=[]
      for x4 in xxx:
        y4=int(x4 * 100)
        yyy.append(y4)
      yy.append(yyy)
    yy=np.array(yy)    
    zz=cv2.fitEllipse(yy)
    z1=zz[1][0]/(2*100)
    z2=zz[1][1]/(2*100)
    ratio = 0.0
    if z1>z2 :
      ratio = z2 / z1
    else:
      ratio = z1 / z2

    fw.write(x+' '+str(ratio)+'\n')

    fr.close()
fw.close()


