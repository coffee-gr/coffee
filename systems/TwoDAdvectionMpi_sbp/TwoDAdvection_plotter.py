#!/usr/bin/env python
# encoding: utf-8


"""
TwoDAdvection_plotter.py

Created by Chris Stevens 2023
"""

# Import Python libraries
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from coffee.actions import Prototype

class Plotter(Prototype):
    def __init__(self, system, frequency = 1, xlim = (-1,1), ylim = (-1,1), \
                 delay = 0.0):
        super(Plotter, self).__init__(frequency)
        self.system = system
        self.delay  = delay
        plt.ion()
        fig = plt.figure(0,figsize=(10,10))
        ax = fig.add_subplot(111)
        
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)  
        ax.grid(True)
        self.axes = ax
        fig.show()
        self.fig = fig

    def _doit(self, it, U):
        x   = U.domain.meshes[0]
        y   = U.domain.meshes[1]
        Psi = U.data[0]
        
        plt.clf()
        self.axes = self.fig.add_subplot(111)
        self.axes.set_title("Iteration: %d, Time: %f" % (it,U.time))
        plt.ioff()
        plt.xlabel('x')
        plt.ylabel('y')

        cplot = self.axes.contourf(x, y, Psi, cmap=plt.cm.bone, \
                                   levels=np.linspace(-1,1,20))
        self.fig.colorbar(cplot)

        plt.ion()
        plt.figure(0)
        plt.draw()
        plt.pause(self.delay)