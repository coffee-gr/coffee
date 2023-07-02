#!/usr/bin/env python
# encoding: utf-8


"""
OneDAdvection_plotter.py

Created by Chris Stevens 2023
"""

# Import Python libraries
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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
        self.line = ax.add_line(Line2D(xlim, ylim, label=r"$\Psi$"))
        
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)       
        ax.grid(True)
        self.axes = ax
        fig.show()

    def _doit(self, it, U):
        x   = U.domain.meshes[0]
        Psi = U.data[0]
        
        self.axes.set_title("Iteration: %d, Time: %f" % (it,U.time))
        plt.ioff()

        self.line.set_xdata(x)
        self.line.set_ydata(Psi)
        self.line.set_color('b')
       
        self.axes.legend(loc = 'lower center',prop={'size': 24},ncol=3)

        plt.ion()
        plt.figure(0)
        plt.draw()
        plt.pause(self.delay)