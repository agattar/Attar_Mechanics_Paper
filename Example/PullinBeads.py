# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 15:38:54 2022

@author: Goktug
"""

# def IndexChanger(beads_pull, coor):
   
beads_pull = 80
coor = 28
import numpy as np
import random

global b_list
global Each_Chr
global c_het
global ind

global bondnumb
global bondarray
global beadname
global bondnumbers


file = open("data.crosslinked","r")
lines = file.readlines()
file.close()

lines[3] = "10 atom types\n"
lines[23] = "9 1\n"
lines.insert(24,"10 1\n\n")
lines.insert(35,"9 1 1\n10 1 1\n\n")

ind = int(lines[4][:5])
atom_num = int(lines[2][:5])

for k,i in enumerate(lines):
    if 'Bonds' in i:
        bondline = k+2
        
bondline = bondline + ind 
for k,i in enumerate(lines):
    if 'Atoms' in i:
        atomline = k+2 


all_list = []
for i in lines[atomline:atomline+atom_num]:
    i = i.split(' ')
    if i[1] != "1":
        all_list.append(i)
        
    
c_list = []
for i in lines[atomline:atomline+atom_num]:
    i = i.split(' ')
    if i[1] == "1":
        c_list.append(i)


    
dnaloc = np.array(c_list);dnaloc = dnaloc[:,:6]  
dnaloc = dnaloc.astype(float)      

dnaloc = dnaloc[dnaloc[:, 0].argsort()]
    
pos_y = np.array(())
neg_y = np.array(())

for k in dnaloc:
    if float(k[4]) > coor:
        pos_y = np.append(k[0],pos_y)        
    if float(k[4]) < -coor:
        neg_y = np.append(k[0],neg_y) 
        
for t in range(beads_pull-1):
    s = random.choice(pos_y)
    t = np.where(pos_y == s)
    for k in dnaloc:
        if s == k[0]:
            k[1] = "10"
            k[2] = "10"
            pos_y = np.delete(pos_y,t)

for t in range(beads_pull-1):
    s = random.choice(neg_y)
    t = np.where(neg_y == s)
    for k in dnaloc:
        if s == k[0]:
            k[1] = "9"
            k[2] = "9"
            neg_y = np.delete(neg_y,t)
final_list = []

for s in dnaloc:
    final_list.append(str(int(s[0])) + " " + str(int(s[1])) + " " +str(int(s[2])) +" "+str((s[3]))+" "+
                      str((s[4]))+ " "+ str((s[5]))+ "\n")

for s in all_list:
    final_list.append(str(int(s[0])) + " " + str(int(s[1])) + " " +str(int(s[2])) +" "+str((s[3]))+" "+
                      str((s[4]))+ " "+ str((s[5]))+ "\n")
    

    
del lines[atomline:atomline+len(final_list)]

lines.insert(atomline,final_list)

fl = open("data.pull",'w+')
     

for line in lines:
    fl.write("".join(line))
    
fl.close()

# if IndexChanger:        
#     IndexChanger(80,28)
#     print("Changed")
# else:
#     pass