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


import logging


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
        self.log = logging.getLogger("IBVP")
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
        self.log.info("Running system %s" % str(self.theSystem))
        self.log.info("Grid = %s" % str(self.theGrid))
        self.log.info("Stepsizes = %s" % repr(u.domain.step_sizes))
        if __debug__:
            self.log.debug("self.actions is %s" % repr(self.theActions))
            self.log.debug("Initial data is = %s" % repr(u))

        advance = self.theSolver.advance
        computation_valid = True
        while computation_valid and t < tstop:
            if __debug__:
                self.log.debug("Beginning new iteration")

            # Check against maxIteration
            if self.iteration > self.maxIteration:
                self.log.warning("Maximum number of iterations exceeded")
                break

            dt = self.theSystem.timestep(u)

            # Check dt for size
            if dt < self.minTimestep:
                self.log.error(
                    "Exiting computation: timestep too small dt = %.15f" % dt
                )
                break

            # Check if dt needs to change in order to hit the next thits value.
            timeleft = tstop - t
            if timeleft < dt:
                dt = timeleft
                if not thits:
                    self.log.warning("Final time step: adjusting to dt = %.15f" % dt)
                    computation_valid = False
                else:
                    self.log.warning("Forcing evaluation at time %f" % tstop)
                    tstop = thits.pop()

            if __debug__:
                self.log.debug("Using timestep dt = %f" % dt)

            # Run the actions.
            self._run_actions(t, u)

            if __debug__:
                self.log.debug("About to advance for iteration = %i" % self.iteration)
            try:
                t, u = advance(t, u, dt)
            except OverflowError as e:
                # OverflowErrors arn't always appropirately handled. So we
                # do it ourselves here.
                print("Overflow error({0}): {1}".format(e.errno, e.strerror))
                computation_valid = False

            # If we're using an mpi enable grid, this ensures that all
            # processes have gotten to the same point before continuing the
            # simulation.
            u.barrier()

            # On to the next iteration.
            self.iteration += 1
            if __debug__:
                self.log.debug("time slice after advance = %s" % repr(u))

        # end (while)

        # Run the actions once more before exiting.
        self._run_actions(t, u)

        # This statement might be unnecessary. In principle it ensures that all
        # processes are about to complete the current simulation before exit
        # occurs.
        u.barrier()
        self.log.info(
            "Finished computation at time %f for iteration %i" % (t, self.iteration)
        )
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
        if __debug__:
            self.log.debug("Running actions")
        tslice = u.collect_data()
        if tslice is not None:
            if __debug__:
                self.log.debug("tslice is not None. Computing if actions will run")
            actions_do_actions = [
                action.will_run(self.iteration, tslice) for action in self.theActions
            ]
            if any(actions_do_actions):
                if __debug__:
                    self.log.debug("Some actions will run")
                for i, action in enumerate(self.theActions):
                    if actions_do_actions[i]:
                        if __debug__:
                            self.log.debug(
                                "Running action %s at iteration %i"
                                % (str(action), self.iteration)
                            )
                        action(self.iteration, tslice)
            if __debug__:
                self.log.debug("All actions done")
