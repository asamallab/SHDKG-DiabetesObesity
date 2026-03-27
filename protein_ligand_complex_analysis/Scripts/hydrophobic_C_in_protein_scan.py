
from math import *
import math
import re
import numpy as np
class atom:
    aid=0    
    atype='' 
    x=0.0  
    y=0.0    
    z=0.0    
    rid=0    
    rtype='' 
    model=[]
    chainid=''
    symb=''

def getlen(atm1,atm2):
    dist=sqrt(pow(atm1.x-atm2.x,2)+pow(atm1.y-atm2.y,2)+pow(atm1.z-atm2.z,2)) 
    return dist

def getangle(atm1,atm2,atm3):
    dist1=sqrt(pow(atm1.x-atm2.x,2)+pow(atm1.y-atm2.y,2)+pow(atm1.z-atm2.z,2)) 
    dist2=sqrt(pow(atm3.x-atm2.x,2)+pow(atm3.y-atm2.y,2)+pow(atm3.z-atm2.z,2)) 
    dotp=(atm1.x-atm2.x)*(atm3.x-atm2.x)+(atm1.y-atm2.y)*(atm3.y-atm2.y)+(atm1.z-atm2.z)*(atm3.z-atm2.z) 
    angle=acos(dotp/(dist1*dist2))*180/pi 
    return angle
def getvanderwaalsum(atm1,atm2):
    if atm1.symb=='N':
        r1=1.55
    if atm1.symb=='O':
        r1=1.52
    if atm1.symb=='S':
        r1=1.8
    if atm1.symb=='C':
        r1=1.7
    if atm1.symb=='Se' or atm1.symb=='SE':
        r1=1.9

    if atm2.symb=='N':
        r2=1.55
    if atm2.symb=='O':
        r2=1.52
    if atm2.symb=='S':
        r2=1.8
    if atm2.symb=='F':
        r2=1.47
    if atm2.symb=='CL' or atm2.symb=='Cl':
        r2=1.75
    if atm2.symb=='BR' or atm2.symb=='Br':
        r2=1.85
    if atm2.symb=='B':
        r2=2
    if atm2.symb=='HG' or atm2.symb=='Hg':
        r2=1.55
    if atm2.symb=='I':
        r2=1.98
    if atm2.symb=='NA' or atm2.symb=='Na':
        r2=2.27
    if atm2.symb=='P':
        r2=1.8
    if atm2.symb=='C':
        r2=1.7
    van=r1+r2+(0.1*(r1+r2))
    return van

def getangledihedral(atm1,atm2,atm3,atm4):
    ab=np.zeros(3)
    bc=np.zeros(3)
    cd=np.zeros(3)
    p=[]
    q=[]
    ab[0]=atm2.x-atm1.x
    ab[1]=atm2.y-atm1.y
    ab[2]=atm2.z-atm1.z
    bc[0]=atm3.x-atm2.x
    bc[1]=atm3.y-atm2.y
    bc[2]=atm3.z-atm2.z
    cd[0]=atm4.x-atm3.x
    cd[1]=atm4.y-atm3.y
    cd[2]=atm4.z-atm3.z
    p.append(ab[1]*bc[2]-ab[2]*bc[1])
    p.append(ab[2]*bc[0]-ab[0]*bc[2])
    p.append(ab[0]*bc[1]-ab[1]*bc[0])
    q.append(bc[1]*cd[2]-bc[2]*cd[1])
    q.append(bc[2]*cd[0]-bc[0]*cd[2])
    q.append(bc[0]*cd[1]-bc[1]*cd[0])


    r1=0
    r2=0
    dp=0
    dpcd=0
    for i in range(0,3):
        r1 += math.pow(p[i],2)
        r2 += math.pow(q[i],2)
        dp += p[i]*q[i]
        dpcd += p[i]*cd[i]

    dih=(dpcd/abs(dpcd))*math.acos(dp/(math.sqrt(r1)*math.sqrt(r2)))*180/math.pi
    

    return dih

filetxt=open('filelist.txt') 
txt_lines=filetxt.read().split('\n') 
filetxt.close()
fileout=open('out_hydrophobic_C_in_protein_scan.txt','w')
f1=open('error_hydrophobic_C_in_protein_scan.txt','w')
intr=[]
lenlines=len(txt_lines)
for ppp in range(lenlines):
    filename=txt_lines[ppp]
    if filename=='':
        continue
    print('%.2f'%((ppp+1)*100.0/(lenlines-1))+'% ('+str(ppp+1)+'/'+str(lenlines-1)+')  Executing for:'+filename)
    fileo=open(filename,'r')
    lines=fileo.read().split('\n')
    fileo.close()
    L=[]
    P=[]
    
    modelno=[]

 
    try:
        for ln in lines:
            if len(ln)>=6 and (ln[0:4]=='ATOM' or ln[0:6]=='HETATM'):
                if ln[0:4]=='ATOM' :
                    atm=atom()
                    atm.aid=int(ln[6:11]) 
                    atm.atype=ln[12:16].strip() 
                    atm.rtype=ln[17:20].strip() 
                    atm.chainid=ln[21]
                    atm.rid=int(ln[22:26]) 
                    atm.x=float(ln[30:38]) 
                    atm.y=float(ln[38:46]) 
                    atm.z=float(ln[46:54]) 
                    atm.model=modelno
                    atm.symb=ln[76:78].strip()
                    #print(atm.aid)
                    
                    if (atm.symb=='C')and (modelno=='1' or modelno=='A' or modelno==[]) :
                        P.append(atm)  
                    
                if ln[0:6]=='HETATM' :
                    atm=atom()
                    atm.aid=int(ln[6:11]) 
                    atm.atype=ln[12:16].strip() 
                    atm.rtype=ln[17:20].strip() 
                    atm.chainid=ln[21]
                    atm.rid=int(ln[22:26]) 
                    atm.x=float(ln[30:38]) 
                    atm.y=float(ln[38:46]) 
                    atm.z=float(ln[46:54]) 
                    atm.model=modelno
                    atm.symb=ln[76:78].strip()
                    #print(atm.aid)     
                    if (atm.symb=='S' or atm.symb=='F' or atm.symb=='Br' or atm.symb=='BR' or atm.symb=='CL' or atm.symb=='Cl' or atm.symb=='I')and (modelno=='1' or modelno=='A' or modelno==[]) :
                        L.append(atm)      
            elif len(ln)>=5 and ln[0:5]=='MODEL':
                modelno=int(ln[12:])

    except:
        f1.write(filename+'\n')

   

    for l in range(len(L)):
        for p in range(len(P)): 
            if getlen(L[l],P[p])<= 4.0:
                  
                intr.append([])
                intr[len(intr)-1].append(filename)                                                     
                intr[len(intr)-1].append(P[p].chainid)
                intr[len(intr)-1].append(P[p].rtype)                
                intr[len(intr)-1].append(P[p].rid)
                intr[len(intr)-1].append(P[p].aid)
                intr[len(intr)-1].append(P[p].atype)
                intr[len(intr)-1].append(L[l].chainid)
                intr[len(intr)-1].append(L[l].rtype)
                intr[len(intr)-1].append(L[l].rid)
                intr[len(intr)-1].append(L[l].aid)
                intr[len(intr)-1].append(L[l].atype)     
                intr[len(intr)-1].append(getlen(P[p],L[l]))           
                

    L=[]
    P=[]
    
    for line in intr:
        for xxd in line:
            fileout.write(str(xxd))
            fileout.write('\t')
        fileout.write('\n')
    intr=[]
    fileout.close()
    fileout=open('out_hydrophobic_C_in_protein_scan.txt','a')
fileout.close()
f1.close()
