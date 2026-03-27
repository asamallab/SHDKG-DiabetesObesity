
from math import *
import math
import re
import numpy as np
import sys
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
    armi=0
    armtype=''

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
def getanglebetweenplanes(atm1,atm2,atm3,atm4,atm5,atm6):
    p1=(atm2.y-atm1.y)*(atm3.z-atm2.z)-(atm2.z-atm1.z)*(atm3.y-atm2.y)
    p2=(atm2.z-atm1.z)*(atm3.x-atm2.x)-(atm2.x-atm1.x)*(atm3.z-atm2.z)
    p3=(atm2.x-atm1.x)*(atm3.y-atm2.y)-(atm2.y-atm1.y)*(atm3.x-atm2.x)
    k1=(atm5.y-atm4.y)*(atm6.z-atm5.z)-(atm5.z-atm4.z)*(atm6.y-atm5.y)
    k2=(atm5.z-atm4.z)*(atm6.x-atm5.x)-(atm5.x-atm4.x)*(atm6.z-atm5.z)
    k3=(atm5.x-atm4.x)*(atm6.y-atm5.y)-(atm5.y-atm4.y)*(atm6.x-atm5.x)
    p=sqrt(pow(p1,2)+pow(p2,2)+pow(p3,2))
    k=sqrt(pow(k1,2)+pow(k2,2)+pow(k3,2))
    angle=(acos((p1*k1+p2*k2+p3*k3)/(p*k)))*180/math.pi
    return angle

def getvanderwaalsum(atm1,atm2):
    if atm1.symb=='N':
        r1=1.55
    if atm1.symb=='O':
        r1=1.52
    if atm1.symb=='S':
        r1=1.8
    if atm1.symb=='Se' or atm1.symb=='SE':
        r1=1.9
    if atm1.symb=='C':
        r1=1.7
    if atm1.symb=='H':
        r1=1.2

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
    if atm2.symb=='H':
        r2=1.2
    van=r1+r2
    return van

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
 
    return False

filetxt=open('filelist.txt') 
txt_lines=filetxt.read().split('\n') 
filetxt.close()
fileout=open('out_aromatic_interaction_f2h_7mem_scan.txt','w')
f1=open('error_aromatic_interaction_f2h_7mem_scan.txt','w')
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
    D=[]
    H=[]
    A=[]
    y=[]
    X1=[]
    X2=[]
    Z=[]
    C1=[]
    C2=[]
    C3=[]
    C4=[]
    C5=[]
    C6=[]
   
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
                    if (atm.rtype=='PHE' or atm.rtype=='TYR')  and (modelno=='1' or modelno=='A' or modelno==[]) :
                        if atm.atype=='CG':
                            C1.append(atm)
                        if atm.atype=='CD1':
                            C2.append(atm)
                        if atm.atype=='CD2':
                            C3.append(atm)
                        if atm.atype=='CE1':
                            C4.append(atm)
                        if atm.atype=='CE2':
                            C5.append(atm)
                        if atm.atype=='CZ':
                            C6.append(atm)
                    if (atm.rtype=='TRP')  and (modelno=='1' or modelno=='A' or modelno==[]) :
                        if atm.atype=='CD2':
                            C1.append(atm)
                        if atm.atype=='CE2':
                            C2.append(atm)
                        if atm.atype=='CE3':
                            C3.append(atm)
                        if atm.atype=='CZ2':
                            C4.append(atm)
                        if atm.atype=='CZ3':
                            C5.append(atm)
                        if atm.atype=='CH2':
                            C6.append(atm)
                    
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
                    if (atm.symb=='C' or atm.symb=='N' or atm.symb=='O' or atm.symb=='S') and (modelno=='1' or modelno=='A' or modelno==[]) :
                        D.append(atm) 
                        #print(filename[6:-4])     
            elif len(ln)>=5 and ln[0:5]=='MODEL':
                modelno=int(ln[12:])

    except:
        f1.write(filename+'\n')
    
    #name=filename[6:-4]+"_docked_pose_with_CathepsinB.aromaticrings.txt"
    #name='AromaticRings/'+filename.strip().split('/')[-1].replace('.pdb','.aromaticrings.txt')
    name='AromaticRings/'+filename.strip().split('/')[-1].replace('_combined.pdb','.aromaticrings.txt')
    f2=open(name,'r+')
    for x in f2:
        y.append(x.split('\t'))
        #print('hi')
    for i in range(0,len(y)):
        if is_number(y[i][0]) == True and float(y[i][2])==7:
            for d in range(len(D)):
                if D[d].x==float(y[i][4]) and D[d].y==float(y[i][5]) and D[d].z==float(y[i][6]):
                    atm=atom()
                    atm.armi=float(y[i][0])
                    atm.x=D[d].x 
                    atm.y=D[d].y 
                    atm.z=D[d].z
                    atm.armtype==(y[i][1])
                    atm.aid=D[d].aid 
                    atm.atype=D[d].atype 
                    atm.rtype=D[d].rtype 
                    atm.rid=D[d].rid
                    atm.chainid=D[d].chainid
                    Z.append(atm)
                    #print(atm.x)

    #print(len(Z))
    
 
    for z1 in range(len(Z)):
        for z2 in range(len(Z)):
            for z3 in range(len(Z)):
                for z4 in range(len(Z)):
                    for z5 in range(len(Z)):
                        for z6 in range(len(Z)):
                            for z7 in range(len(Z)):
                                if Z[z1].armi==Z[z2].armi==Z[z3].armi==Z[z4].armi==Z[z5].armi==Z[z6].armi==Z[z7].armi and Z[z1].aid>Z[z2].aid>Z[z3].aid>Z[z4].aid>Z[z5].aid>Z[z6].aid>Z[z7].aid:
                                    atm=atom()
                                    atm.armi=Z[z1].armi
                                    atm.x=(Z[z1].x+Z[z2].x+Z[z3].x+Z[z4].x+Z[z5].x+Z[z6].x+Z[z7].x)/7
                                    atm.y=(Z[z1].y+Z[z2].y+Z[z3].y+Z[z4].y+Z[z5].y+Z[z6].y+Z[z7].y)/7
                                    atm.z=(Z[z1].z+Z[z2].z+Z[z3].z+Z[z4].z+Z[z5].z+Z[z6].z+Z[z7].z)/7
                                    atm.armtype==Z[z1].armtype
                                    #atm.aid=Z[z1].aid 
                                    #atm.atype=Z[z1].atype 
                                    atm.rtype=Z[z1].rtype 
                                    atm.rid=Z[z1].rid
                                    atm.chainid=Z[z1].chainid
                                    X2.append(atm)
                                
                                   
    if len(C1)==len(C2)==len(C3)==len(C4)==len(C5)==len(C6):
        for c1 in range(len(C1)):
            atm=atom()
            #atm.aid=int(ln[6:11]) 
            #atm.atype=ln[12:16].strip() 
            atm.rtype=C1[c1].rtype 
            atm.chainid=C1[c1].chainid
            atm.rid=C1[c1].rid
            atm.x=(C1[c1].x+C2[c1].x+C3[c1].x+C4[c1].x+C5[c1].x+C6[c1].x)/6  
            atm.y=(C1[c1].y+C2[c1].y+C3[c1].y+C4[c1].y+C5[c1].y+C6[c1].y)/6 
            atm.z=(C1[c1].z+C2[c1].z+C3[c1].z+C4[c1].z+C5[c1].z+C6[c1].z)/6 
            X1.append(atm)

    for c1 in range(len(C1)):
        for c2 in range(len(C2)):
            for x1 in range(len(X1)):
                if C1[c1].rid==C2[c2].rid==X1[x1].rid and C1[c1].chainid==C2[c2].chainid==X1[x1].chainid:
                    for z1 in range(len(Z)):
                        for z2 in range(len(Z)):
                            for x2 in range(len(X2)):
                                if Z[z1].rid==Z[z2].rid==X2[x2].rid and Z[z1].chainid==Z[z2].chainid==X2[x2].chainid and Z[z1].aid>Z[z2].aid  and getlen(X1[x1],X2[x2])<=5.5 and getanglebetweenplanes(C1[c1],C2[c2],X1[x1],Z[z1],Z[z2],X2[x2])<=120 and getanglebetweenplanes(C1[c1],C2[c2],X1[x1],Z[z1],Z[z2],X2[x2])>=60:
                                    intr.append([])
                                    intr[len(intr)-1].append(filename)                                                     
                                    intr[len(intr)-1].append(X1[x1].chainid)
                                    intr[len(intr)-1].append(X1[x1].rtype)                
                                    intr[len(intr)-1].append(C1[c1].rid)
                                    intr[len(intr)-1].append(X2[x2].chainid)
                                    intr[len(intr)-1].append(X2[x2].rtype)                
                                    intr[len(intr)-1].append(Z[z2].rid)
                                    intr[len(intr)-1].append(X2[x2].armi)
                                    intr[len(intr)-1].append(getlen(X1[x1],X2[x2]))           
                                    intr[len(intr)-1].append(getanglebetweenplanes(C1[c1],C2[c2],X1[x1],Z[z1],Z[z2],X2[x2]))

                            
                


              
 

    D=[]
    H=[]
    A=[]
    C1=[]
    C2=[]
    C3=[]
    C4=[]
    C5=[]
    C6=[]
    X1=[]
    X2=[]
    Z=[]
    for line in intr:
        for xxd in line:
            fileout.write(str(xxd))
            fileout.write('\t')
        fileout.write('\n')
    intr=[]
    fileout.close()
    fileout=open('out_aromatic_interaction_f2h_7mem_scan.txt','a')
fileout.close()
f1.close()
