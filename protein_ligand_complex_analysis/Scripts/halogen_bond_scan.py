
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

def getcxlength(atm1,atm2):
    if atm1.symb=='C' and atm2.symb=='F' :
        r1=1.45
    if atm1.symb=='C' and (atm2.symb=='Cl' or atm2.symb=='CL'):
        r1=1.85
    if atm1.symb=='C' and (atm2.symb=='Br' or atm2.symb=='BR'):
        r1=2.04
    if atm1.symb=='C' and (atm2.symb=='I'):
        r1=2.24
     
    return r1


filetxt=open('filelist.txt') 
txt_lines=filetxt.read().split('\n') 
filetxt.close()
fileout=open('out_halogen_bond_scan.txt','w')
f1=open('error_halogen_bond_scan.txt','w')
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
    C=[]
    X=[]
    Y=[]
    B=[]
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
                    
                    if (atm.atype=='N' or atm.atype=='O'  ) and (modelno=='1' or modelno=='A' or modelno==[]) :
                        Y.append(atm) 
                    if (atm.atype=='C' or atm.atype=='CA'  ) and (modelno=='1' or modelno=='A' or modelno==[]) :
                        B.append(atm)  
                    if (atm.rtype=='ASN' or atm.rtype=='ASP' ) and (modelno=='1' or modelno=='A' or modelno==[]) :
                        if atm.atype=='OD1' or atm.atype=='ND2':
                            Y.append(atm)  
                        if atm.atype=='CG':
                            B.append(atm) 
                    if (atm.rtype=='GLN' or atm.rtype=='GLU' ) and (modelno=='1' or modelno=='A' or modelno==[]) :
                        if atm.atype=='OE1' or atm.atype=='NE2':
                            Y.append(atm)  
                        if atm.atype=='CD':
                            B.append(atm)  
                    if (atm.rtype=='ARG' ) and (modelno=='1' or modelno=='A' or modelno==[]) :
                        if atm.atype=='NE' or atm.atype=='NH1' or atm.atype=='NH2':
                            Y.append(atm)  
                        if atm.atype=='CZ':
                            B.append(atm)
                    if (atm.rtype=='HIS' ) and (modelno=='1' or modelno=='A' or modelno==[]) :
                        if atm.atype=='ND1' or atm.atype=='NE2' :
                            Y.append(atm)  
                        if atm.atype=='CE1':
                            B.append(atm)  
                    if (atm.rtype=='MET' ) and (modelno=='1' or modelno=='A' or modelno==[]) :
                        if atm.atype=='SD' :
                            Y.append(atm)  
                        if atm.atype=='CE':
                            B.append(atm)  
                    if (atm.rtype=='SER' ) and (modelno=='1' or modelno=='A' or modelno==[]) :
                        if atm.atype=='OG' :
                            Y.append(atm)  
                        if atm.atype=='CB':
                            B.append(atm) 
                    if (atm.rtype=='THR' ) and (modelno=='1' or modelno=='A' or modelno==[]) :
                        if atm.atype=='OG1' :
                            Y.append(atm)  
                        if atm.atype=='CB':
                            B.append(atm)
                    if (atm.rtype=='TRP' ) and (modelno=='1' or modelno=='A' or modelno==[]) :
                        if atm.atype=='NE1' :
                            Y.append(atm)  
                        if atm.atype=='CD1':
                            B.append(atm)
                    if (atm.rtype=='TYR' ) and (modelno=='1' or modelno=='A' or modelno==[]) :
                        if atm.atype=='OH' :
                            Y.append(atm)  
                        if atm.atype=='CZ':
                            B.append(atm)  
                    if (atm.rtype=='CYS' ) and (modelno=='1' or modelno=='A' or modelno==[]) :
                        if atm.atype=='SG' :
                            Y.append(atm)  
                        if atm.atype=='CB' :
                            B.append(atm)
                        
                    
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
                    if (atm.symb=='F' or atm.symb=='Cl' or atm.symb=='CL' or atm.symb=='Br' or atm.symb=='BR' or atm.symb=='I' ) and (modelno=='1' or modelno=='A' or modelno==[]) :
                        X.append(atm) 
                    if (atm.symb=='C' ) and (modelno=='1' or modelno=='A' or modelno==[]) :
                        C.append(atm)      
            elif len(ln)>=5 and ln[0:5]=='MODEL':
                modelno=int(ln[12:])

    except:
        f1.write(filename+'\n')

   

    for y in range(len(Y)):
        for b in range(len(B)):
            if B[b].rid==Y[y].rid and B[b].chainid==Y[y].chainid and getlen(B[b],Y[y])<= 1.9:
                for c in range(len(C)):
                    for x in range(len(X)):
                        if C[c].rid==X[x].rid and C[c].chainid==X[x].chainid and getlen(C[c],X[x])<=getcxlength(C[c],X[x]) and getlen(X[x],Y[y])<= 3.7 and getangle(C[c],X[x],Y[y])>=90 and getangle(X[x],Y[y],B[b])>=90 :
                    
                            intr.append([])
                            intr[len(intr)-1].append(filename)                                                     
                            intr[len(intr)-1].append(Y[y].chainid)
                            intr[len(intr)-1].append(Y[y].rtype)                
                            intr[len(intr)-1].append(Y[y].rid)
                            intr[len(intr)-1].append(Y[y].aid)
                            intr[len(intr)-1].append(Y[y].atype)
                            intr[len(intr)-1].append(B[b].atype)
                            intr[len(intr)-1].append(B[b].aid)
                            intr[len(intr)-1].append(C[c].chainid)
                            intr[len(intr)-1].append(C[c].rtype)
                            intr[len(intr)-1].append(C[c].rid)
                            intr[len(intr)-1].append(C[c].aid)
                            intr[len(intr)-1].append(C[c].atype)
                            intr[len(intr)-1].append(X[x].atype)  
                            intr[len(intr)-1].append(X[x].aid)   
                            intr[len(intr)-1].append(getlen(X[x],Y[y]))           
                            intr[len(intr)-1].append(getangle(X[x],Y[y],B[b]))
                            intr[len(intr)-1].append(getangle(C[c],X[x],Y[y]))

    C=[]
    X=[]
    Y=[]
    B=[]
    for line in intr:
        for xxd in line:
            fileout.write(str(xxd))
            fileout.write('\t')
        fileout.write('\n')
    intr=[]
    fileout.close()
    fileout=open('out_halogen_bond_scan.txt','a')
fileout.close()
f1.close()
