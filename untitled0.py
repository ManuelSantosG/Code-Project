#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 11:08:29 2018

@author: ms5717
"""


def SeqN(a):
    l=[]
    n=len(a)
    c=1
    i=0
    k=0
    p=0
    while i < n-1:
        if a[i]==a[i+1]:
            Counter=True;
            for j in range(i,n-1):
                if a[i]==a[j+1]:
                    c+=1
                    k=i+c
                else:
                    k=i+c
            l.append(c)
            l.append(a[i])
            p+=c
            
        else:
            Counter=False
            l.append(1)
            l.append(a[i])
            p+=1

        c=1
        if Counter:
            i+=c
        else:
            i+=1
    if p<n:
        l.append(1)
        l.append(a[n-1])
    return(l,len(l))
SeqN([1])
        
    