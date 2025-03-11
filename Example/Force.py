# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 11:10:04 2022

@author: Asus
"""

import matplotlib.pyplot as plt
import numpy as np
import os


files = ["force.txt"]

zz = np.zeros((201,len(files)))

for k,t in enumerate(files):
    file = open(t)
    file = file.readlines()
    for s,i in enumerate(file[2:]):
        i = i.split(" ")
        zz[s,k] = abs(float(i[2]))
        
np.savetxt("Force_y.txt",zz,delimiter=',')