#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 09:19:44 2018

@author: blake
"""

import math
from scipy import optimize


def fun(x):
    return [x[0]**2-4,
            math.exp(x[1])-math.exp(-3),
            (x[0]-2)**2+(x[1]+3)**2-1e-8 ]


sol = optimize.root(fun, [0, 0], method='lm')
print(sol.x)





