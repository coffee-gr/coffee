"""A module containing actions that interface with GnuPlot."""


from builtins import str
from builtins import range
import time
import os
from coffee.backend import backend as be

# Gnuplot currently only works with Python 2!
try:
    import Gnuplot
except ImportError:
    import sys

    print(sys.exc_info())
    print("Gnuplot currently only works with Python 2!")
    sys.exit()

from coffee.actions import Prototype


class Plotter1D(Prototype):
    """This class provides convient access to GnuPlot's functions for one
    dimensional data.

    It is not
    suitable for animations since the GnuPlot object is persistent between
    calls to this class.

    The simulation itself does not need to be one dimensional. The data to be
    plotted is 'pre-processed' by a function that is passed by the user. In
    this way the user can, for example, specify if the modulus, or argument
    of complex data is to be plotted.
    """

    def __init__(self, system, *args, **kwds):
        """The constructor for Plotter1D objects.

        Plotter1D subclasses actions.Prototype.

        All positional arguments, excluding the
        first (``system``), are passed to the GnuPlot object. These
        arguments should be strings respecting the documentation of GnuPlot.
        They are passed to the Gnuplot object g as g(pos_argument)

        Additional parameters valid for actions.Prototype are valid.

        Parameters
        ==========
        system : system
            The system object.
        args : list of strings
            Passed to gnuplot as arguments.
        delay : float, Optional
            The delay in seconds between plotting events.
        title : string, Optional
            The title of the plot. The default is the time of the plot.
        data_function : function(it, u, system) -> (be.ndarray, be.ndarray), Optional
              A function that takes the interation number,
              a timeslice and the system. It should return a tuple consisting
              of the domain of the data and a two dimensional array. The first
              dimension of the array is considered to be the components. The
              second dimension is considered to give the domain dependence.
        """

        if "start" in kwds:
            start = kwds.pop("start")
        else:
            start = -float("Infinity")
        if "frequency" in kwds:
            frequency = kwds.pop("frequency")
        else:
            frequency = 1
        if "delay" in kwds:
            self.delay = kwds.pop("delay")
        else:
            self.delay = 0.0
        self.title = r"Time %f"
        if "title" in kwds:
            self.title = kwds.pop("title")
        self.system = system
        try:
            if kwds["data_function"] is not None:
                self.datafunc = kwds.pop("data_function")
            else:
                self.datafunc = lambda y, x, z: (x.domain.axes[0], x.data)
        except:
            self.datafunc = lambda y, x, z: (x.domain.axes[0], x.data)
        super(Plotter1D, self).__init__(frequency=frequency, start=start, **kwds)
        self.Device = Gnuplot.Gnuplot()
        g = self.Device
        for arg in args:
            g(arg)

    def _doit(self, it, u):
        g = self.Device
        x, f = self.datafunc(it, u, self.system)
        f = be.atleast_2d(f)
        graphs = []
        g("set title '" + self.title % u.time + "'")
        g.plot(*graphs)
        time.sleep(self.delay)

    def __del__(self):
        del self.Device


class Plotter2D(Prototype):
    """This class provides convient access to GnuPlot's functions for two
    dimensional data.

    It is not
    suitable for animations since the GnuPlot object is persistent between
    calls to this class.

    The simulation itself does not need to be two dimensional. The data to be
    plotted is 'pre-processed' by an function that is passed by the user.
    """

    def __init__(self, *args, **kwds):
        """The constructor for Plotter2D objects.

        Plotter2D subclasses actions.Prototype. Any argument valid
        in the initialiser of actions.Prototype is respected here.

        All positional arguments
        are passed to the GnuPlot object. These
        arguments should be strings respecting the documentation of GnuPlot.
        They are passed to the Gnuplot object g as g(pos_argument)

        The data_function
        should return a tuple consisting
        of the domain of the data and a three dimensional array. The first
        dimension of the array is considered to be a number of components
        each of which is to be plotted. The
        second and third dimensions is considered to give the domain dependence.

        Parameters
        ==========
        system : system
            The system object.
        args : list of strings
            Passed to gnuplot as arguments.
        delay : float, Optional
            The delay in seconds between plotting events.
        components : list of strings, Optional
            The strings are used to label the components being plotted.
        title : string, Optional
            The title of the plot. The default is the time of the plot.
        data_function : function(it, u, system) -> (be.ndarray, be.ndarray), Optional
              A function that takes the interation number,
              a timeslice and the system. It should return a tuple consisting
              of the domain of the data and a two dimensional array. The first
              dimension of the array is considered to be the components. The
              second dimension is considered to give the domain dependence.
        """
        if "delay" in kwds:
            self.delay = kwds.pop("delay")
        else:
            self.delay = 0.0

        if "title" in kwds:
            self.title = kwds.pop("title")
        else:
            self.title = "Iteration %d, Time %f"

        if "components" in kwds:
            self.components = kwds.pop("components")
        else:
            self.components = None

        self.system = kwds.pop("system")
        try:
            if kwds["data_function"] is not None:
                self.datafunc = kwds.pop("data_function")
            else:
                self.datafunc = lambda y, x, z: (x.domain.axes, x.data)
        except:
            self.datafunc = lambda y, x, z: (x.domain.axes, x.data)

        self.device = Gnuplot.Gnuplot()
        for arg in args:
            self.device(arg)
        super(Plotter2D, self).__init__(**kwds)

    def _doit(self, it, u):
        x, f = self.datafunc(it, u, self.system)
        f = be.atleast_2d(f)
        graphs = []
        self.device('set title "%s" enhanced' % self.title % (it, u.time))
        for i in range(f.shape[0]):
            if self.components is None:
                graphs += [
                    Gnuplot.GridData(
                        f[i, :, :],
                        xvals=x[0],
                        yvals=x[1],
                        filename="deleteme.gp",
                        title="Component %i" % i,
                        binary=0,
                    )
                ]
            else:
                graphs += [
                    Gnuplot.GridData(
                        f[i, :, :],
                        xvals=x[0],
                        yvals=x[1],
                        filename="deleteme.gp",
                        title=self.components[i],
                        binary=0,
                    )
                ]
        self.device.splot(*graphs)
        time.sleep(self.delay)

    def __del__(self):
        try:
            os.remove("deleteme.gp")
        except:
            print(
                "The file deleteme.gp was unable to be deleted. Please do so manually."
            )
        del self.device
