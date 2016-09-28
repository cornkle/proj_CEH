# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 17:27:16 2016

@author: cornkle
"""

import os, fnmatch
import math

# Recursively locates PATTERN under ROOT_PATH for finding the absolute path to files
def rlocate(pattern, root_path):
    strg=[]
    for path, dirs, files in os.walk(os.path.abspath(root_path)):
        for filename in fnmatch.filter(files, pattern):
            strg.append(os.path.join(path, filename))
    return strg
    
    
#def locate(pattern, root_path):
#    strg=[]
#    for file in os.listdir(root_path):
#        if file.endswith('.txt'):
#            strg.append(os.path.join(root_path, file))
#    return strg    
    
def locate(pattern, root_path):
    strg=[]
    llist=os.listdir(root_path)
    llist.sort()
    for file in llist:
        #print(file)
        if file.endswith(pattern):
            strg.append(os.path.join(root_path, file))
    return strg    

def distance(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2)
            
            
            


