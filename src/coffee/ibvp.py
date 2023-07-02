#!/usr/bin/env python
# encoding: utf-8
"""
This is the class that handles the iterative step and actions of the
simulation. 

The name ibvp comes from the original design of this class as the
iterative step in an initial boundary value problem solver. However, this class
is general enough to handle numerical techniques that have, at their highest
level, some iterative method.
"""

class IBVP:
    """Handles computation of the solution to initial boundary value problems.

    This class manages the interative process of simulations and interleaves
    interative steps with calls to a list of actions. Information about events
    during the simulation are writen to logging.getLogger("IBVP").
    """

    theActions = None
    theGrid = None
    theSolver = None
    iteration = 0
    maxIteration = None

    def __init__(self, sol, eqn, grid, action=[], maxIteration=10000, minTimestep=1e-8):
        """IBVP constructor.

        Parameters
        ----------

        sol : solver.Solver
              An ODE solver.

        eqn : system.System
              The system of equations whose solution will be computed.

        grid : grid.Grid
               The spatial grid over which the values of the solution to the system
               will be computed.

        action : list, actions.Prototype
                 Elements of the list must be instances of actions.Prototype. Each
                 action will be run after initialisation of the timeslice and after
                 each new timeslice is calculated.

        maxIteraction : int, Optional
                 Specifies the maximum number of timesteps.

        minTimestep : float, Optinoal
                 Specifies the minimum time step.
        """
        self.theSolver = sol
        self.theSystem = eqn
        self.maxIteration = maxIteration
        self.theGrid = grid
        self.theActions = action
        self.minTimestep = minTimestep

    def run(self, tstart, tstop=float("inf"), thits=None):
        """Go for it! Starts the simulation.

        Parameters
        ----------

        tstart : float
            A number giving the time at the initial interation.

        tstop : float
            A number giving the time to finish the simulation at.

        thits : list, float
            A list of times that the simulation is forced to hit exactly.

        """
        # Set up thits
        if thits is None:
            thits = []
        if tstop not in thits:
            thits += [tstop]

        # Order the list of times, and ensure that they are popped from smallest
        # to largest.
        thits = sorted(thits)
        thits.reverse()
        tstop = thits.pop()

        # Set start time.
        t = tstart

        # Get initial data and configure timeslices for multiple processors.
        u = self.theSystem.initial_data(t, self.theGrid)

        advance = self.theSolver.advance
        computation_valid = True
        while computation_valid and t < tstop:
            # Check against maxIteration
            if self.iteration > self.maxIteration:
                break

            dt = self.theSystem.timestep(u)

            # Check dt for size
            if dt < self.minTimestep:
                break

            # Check if dt needs to change in order to hit the next thits value.
            timeleft = tstop - t
            if timeleft < dt:
                dt = timeleft
                if not thits:
                    computation_valid = False
                else:
                    tstop = thits.pop()

            # Run the actions.
            self._run_actions(t, u)

            try:
                t, u = advance(t, u, dt)
            except OverflowError as e:
                # OverflowErrors arn't always appropirately handled. So we
                # do it ourselves here.
                print("Overflow error({0}): {1}".format(e.errno, e.strerror))
                computation_valid = False

            # Dirichlet boundary conditions
            try:
                u = self.theSystem.dirichlet_boundary(u)
            except AttributeError:
                pass

            # If we're using an mpi enable grid, this ensures that all
            # processes have gotten to the same point before continuing the
            # simulation.
            u.barrier()

            # On to the next iteration.
            self.iteration += 1

        # end (while)

        # Run the actions once more before exiting.
        self._run_actions(t, u)

        # This statement might be unnecessary. In principle it ensures that all
        # processes are about to complete the current simulation before exit
        # occurs.
        u.barrier()

        return u

    def _run_actions(self, t, u):
        """A utility method that ensures that actions are only run after all data is collected across MPI nodes and only if an action is due to run.

        Parameters
        ----------

        t : float
            The current time.

        u: tslices.TimeSlice
            The current timeslice.
        """

        # Ideally u.collect_data() should only be executed if there
        # actions that will run. Because of single process access to
        # actions this causes an issue.
        # Some thought is required to fix this problem.
        tslice = u.collect_data()
        if tslice is not None:
            actions_do_actions = [
                action.will_run(self.iteration, tslice) for action in self.theActions
            ]
            if any(actions_do_actions):
                for i, action in enumerate(self.theActions):
                    if actions_do_actions[i]:
                        action(self.iteration, tslice)