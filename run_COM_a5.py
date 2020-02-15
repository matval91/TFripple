#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 13:14:21 2020

@author: vallar
"""
import COM_a5_class as c5
import sys

if len(sys.argv) == 4:
    fname_a5=sys.argv[1]
    run = sys.argv[2]
    E=float(sys.argv[3])
else:
    fname_a5='/home/vallar/WORK/ASCOT/runs/SA_003/nnb_ripple/production/ascot.h5'
    run='run_1893662328'
    E=500
print('Read input', fname_a5, run, str(E))
debug=False
plot=True
eq=c5.COM_a5(fname_a5, run, E, debug=debug, plot=plot)