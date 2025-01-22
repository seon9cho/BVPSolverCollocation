#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 13:21:36 2018

@author: blake
"""

import numpy as np


def fun(x, y):
    return np.vstack((y[1], -np.exp(y[0])))

def bc(ya, yb):
    return np.array([ya[0], yb[0]])

x = np.linspace(0, 1, 5)

y_guess = np.zeros((2, x.size))

#y_b = np.zeros((2, x.size))
#y_b[0] = 3

from scipy.integrate import solve_bvp
res_a = solve_bvp(fun, bc, x, y_guess)
#res_b = solve_bvp(fun, bc, x, y_b)



# PLOTTING


x_plot = np.linspace(0, 1, 100)
y_plot_a = res_a.sol(x_plot)[0]
#y_plot_b = res_b.sol(x_plot)[0]

import matplotlib.pyplot as plt
#plt.plot(x_plot, y_plot_a, label='y_a')
#plt.plot(x_plot,res_a.sol(x_plot)[1])
#plt.plot(x_plot, y_plot_b, label='y_b')


plt.plot(y_plot_a,res_a.sol(x_plot)[1])

plt.legend()
plt.xlabel("y_1")
plt.ylabel("y_2")
plt.show()









