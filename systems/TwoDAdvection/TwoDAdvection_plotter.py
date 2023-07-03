import numpy as np
from pylab import *

from coffee.actions import Prototype


class Plotter(Prototype):
    """docstring for Plotter"""

    def __init__(
        self, system, frequency=1, xlim=(-1, 1), ylim=(-1, 1), delay=0.0, fignum=0
    ):
        super(Plotter, self).__init__(frequency)

        self.delay = delay
        self.colors = ("b", "g", "r", "c", "m", "y", "k", "orange")
        from numpy import asarray

        ion()
        self.fignum = fignum
        fig = figure(fignum, figsize=(10, 10))
        ax = fig.add_subplot(111)
        self.lines = []

        self.lines.append(ax.add_line(Line2D(xlim, ylim, label=r"$\Psi$")))

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.grid(True)

        self.axes = ax

        self.system = system

        fig.show()

    def _doit(self, it, u):
        self.axes.set_title("Iteration: %d, Time: %f" % (it, u.time))

        ioff()

        x = u.domain.meshes[0][:, 0]
        Psi = u.data[0][:, 20]

        self.lines[0].set_xdata(x)
        self.lines[0].set_ydata(Psi)
        self.lines[0].set_color(self.colors[5])

        self.axes.legend(loc="upper right")

        ion()
        plt.figure(self.fignum)
        draw()
        plt.pause(self.delay)
