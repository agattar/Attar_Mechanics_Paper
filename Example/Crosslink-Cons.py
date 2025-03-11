# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 15:38:54 2022

@author: Goktug
"""

import numpy as np
import random


def Crosslinker(chnum, prob,minbond,maxbond):
    
    file = open("data.crosslinked","r")
    lines = file.readlines()
    file.close()
    
    b_list = np.zeros((chnum*24001+10000,2))
    count = 0
    
    
      
    ind = (lines[4])
    bondind = ""
    for k in ind:
        if k.isdigit():
            bondind += k
            
            
    atom_num = int(lines[2][:5])
    
    for k,i in enumerate(lines):
        if 'Bonds' in i:
            bondline = k+2
            
    bondline = bondline + int(bondind)
    bondind = int(bondind)
    for k,i in enumerate(lines):
        if 'Atoms' in i:
            atomline = k+2 
    
    h_list = []
    for i in lines[atomline:atomline+atom_num]:
        i = i.split(' ')
        if i[1] == "4" or i[1] == "6":
            h_list.append(i)

            
            
    dnaloc = np.array(h_list);dnaloc_f = dnaloc[:,:6].astype(float)        
    dnaloc_f = dnaloc_f[dnaloc_f[:, 0].argsort()]
   
    
    bondnumb = np.random.randint(minbond,maxbond,size=(len(dnaloc),1))
    bondnumbers = np.linspace(1,len(dnaloc),len(dnaloc))
    bondnumbers = np.reshape(bondnumbers,(len(dnaloc),1))
    bondnumbers = np.hstack((bondnumbers,bondnumb))
    
    
    b_index = 0
    bnum = []
    beadname = dnaloc[dnaloc[:, 0].argsort()]
    
    for n in range(len(dnaloc[:,:])):
        
        r = random.random()
        if r <= prob:
            neighbor_list = []
            Extracted = np.subtract(dnaloc_f[n,3:],dnaloc_f[:,3:])
            Multiplied = np.multiply(Extracted,Extracted)
            Multiplied = np.sum(Multiplied,axis=1)  
            Distances = np.sqrt(Multiplied)
        
            Result = np.where((Distances[:] <= 2.50) & (Distances[:] >= 1.00))
    
            for k in Result[0]:
                
                selected_num = int(bondnumbers[k][1])
                
                if ([int(beadname[k][0]),int(beadname[n][0])]==b_list[:,None]).all(-1).any() == False and selected_num > 0:
                    diff = k - n
                    if (diff == 1 or diff == -1):
                        pass
                    else:
                        neighbor_list.append(k)
                else:
                    pass

            while (int(bondnumbers[n,1])) > 0:   
                
                if len(neighbor_list) != 0:
                    chosen = random.choice(neighbor_list)
                    if ([int(beadname[n][0]),int(dnaloc_f[chosen][0])]==b_list[:,None]).all(-1).any() == False and ([int(dnaloc_f[chosen][0]),int(beadname[n][0])]==b_list[:,None]).all(-1).any() == False:
                        bondnumbers[int(n),1] -= 1
                        bondind += 1
                        bnum.append(str(bondind)+" "+"4"+" "+str(int(beadname[n][0]))+" "+str(int(dnaloc_f[chosen][0]))+"\n")
                        b_list[count+b_index,0] = str(beadname[n][0]);b_list[count+b_index,1] = str(dnaloc_f[chosen][0])
                        neighbor_list.remove(chosen)
                        b_index +=1
                    else:
                        neighbor_list.remove(chosen)
                                    
                else:
                    break
                          
        else:
            pass
            
                    
    mean_bondnum = np.subtract(bondnumb,bondnumbers[:,1])   
    mean_bondnum = np.mean(mean_bondnum)
    print(mean_bondnum,"Crosslink Mean Bond Number")
    
    indstr = str(str(bondind)+" "+"bonds"+"\n")
    lines[4] = indstr
    
    bnum.append("\n")
    lines.insert(bondline,bnum)
    
    fl = open("data.crosslinked",'w+')
         
    for line in lines:
        fl.write("".join(line))
        
    fl.close()

if Crosslinker:        
    Crosslinker(8, 0.2, 1, 3)
    print("Crosslinked")
else:
    pass
