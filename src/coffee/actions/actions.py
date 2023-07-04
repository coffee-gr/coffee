from builtins import object
from coffee.settings import be


class Prototype(object):
    """The prototype of all actions.

    This class provides basic functionality that all actions require.

    To subclass this class:
    1) define the method _doit in the subclass
    2) call the Prototypes constructor to ensure that frequency, start and stop
    are properly initialised.
    """

    def __init__(
        self,
        frequency=1,
        start=-float("infinity"),
        stop=float("infinity"),
        thits=None,
        thits_toll=1e-14,
    ):
        """The initialiser for the Prototype action class.

        Parameters
        ==========
        frequency : int, Optional
            The action will only run after the given number of iterations.
        stop : float, Optional
            After this time the action will not run.
        start : float, Optional
            Before this time the action will not run.
        thits : list of floats, Optional
            Times at which the action is guaranteed to run
        thits_toll : float, Optional
            The tollerance for comparison against times in thits.
        """
        self.frequency = frequency
        self.stop = stop
        self.start = start
        if thits is not None:
            self.thits = be.asarray(thits)
        else:
            self.thits = None
        self.thits_toll = thits_toll

    def will_run(self, it, u):
        """Returns true if the action will run for the given iteration and
        time slice.

        Parameters
        ==========
        it : int
            The number of iterations in this simulation.
        u : tslice.TimeSlice
            The current timeslice.

        Returns
        =======
        boolean :
            True if the action will run.
        """
        test = (it % self.frequency) == 0 and self.start <= u.time <= self.stop
        if self.thits is not None:
            test = test and (be.absolute(u.time - self.thits) < self.thits_toll).any()
        return test

    def __call__(self, it, u):
        """If the action will run, runs the action.

        Parameters
        ==========
        it : int
            The number of iterations in this simulation.
        u : tslice.TimeSlice
            The current timeslice.
        """
        if self.will_run(it, u):
            self._doit(it, u)

    def _doit(self, it, u):
        pass


class BlowupCutoff(Prototype):
    """An action that stops the simulation if some component of the computed
    solution is above the given value.
    """

    def __init__(self, cutoff=10, **kwds):
        """The initialiser for the BlowupCutoff action.

        The parameters listed below are only the parameters specific to this
        action. All parameters as given for the initialiser for the
        actions. Prototype actions are also valid.

        Parameters
        ==========
        cutoff: float
            The value that is used to test for ending the simulation.
        """
        super(BlowupCutoff, self).__init__(**kwds)
        self.cutoff = abs(cutoff)

    def above_cutoff(self, u):
        """Returns true if the absolute value of any component of
        the data in the timeslice u is above the absolute value of
        the given cutoff.

        Parameters
        ==========
        u : tslice.TimeSlice
            The timeslice to be checked.

        Returns
        =======
        boolean :
            True if there is a component greater than the cutoff.
        """
        for component in u.fields:
            if be.any(be.abs(component) >= self.cutoff):
                return True
        return False

    def _doit(self, it, u):
        if self.above_cutoff(u):
            raise Exception("Function values are above the cutoff")


class Info(Prototype):
    """Logs the current iterations and data slice.

    Included more as an example than anything else.
    """

    def __init__(self, *args, **kwds):
        """The initialiser for the Info action.

        All parameters as given for the initialiser for the
        actions.Prototype action are valid.
        """
        super(Info, self).__init__(*args, **kwds)

    def _doit(self, it, u):
        pass
