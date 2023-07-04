from builtins import range
import time

# from pylab import *
import pylab
from coffee.settings import be

from coffee.actions import Prototype


class Plotter(Prototype):
    """An action that plots data using matplotlib."""

    def __init__(self, frequency=1, xlim=(-1, 1), ylim=(-1, 1), findex=None, delay=0.0):
        """Initialiser for Plotter.

        Parameters
        ==========
        frequency : int
            The number of iterations per execution of this action
        xlim : a two tuple of floats
            The limits of the x axes.
        ylim : a two tuple of floats
            The limits of the y axes.
        findex : int
            The number of functions to plot
        delay : float
            The delay, in seconds, between repeated executions of this function.
        """
        super(Plotter, self).__init__(frequency)
        if findex is not None:
            self.index = findex
            self.delay = delay
            self.colors = ("b", "g", "r", "c", "m", "y", "k", "coral")
            from numpy import asarray

            pylab.ion()
            fig = figure(1)
            ax = fig.add_subplot(111)
            self.lines = [ax.add_line(Line2D(xlim, ylim)) for k in findex]
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.grid(True)
            self.axes = ax

    def _doit(self, it, u):
        x = u.x
        index = be.asarray(self.index)
        f = u.fields[index]
        mx = be.max(f.flat)
        mn = be.min(f.flat)
        self.axes.set_title("Iteration: %d, Time: %f" % (it, u.time))
        # self.axes.set_xlim(x[0],x[-1])
        # self.axes.set_ylim(mn,mx,auto=True)

        l = len(self.index)
        pylab.ioff()
        for k in range(l):
            line = self.lines[k]
            line.set_xdata(x)
            line.set_ydata(f[k])
            line.set_color(self.colors[k])
        pylab.ion()
        pylab.draw()
        time.sleep(self.delay)
