import logging
from unittest.mock import MagicMock

import pytest

from coffee.actions import actions
from coffee.backend import backend as be

backends = ["numpy"]
# When other backends are added, add them to this list
# e.g. backends = ["numpy", "cupy", "jax"]


@pytest.mark.parametrize("backend", backends)
class TestPrototype:
    def test_will_run(self, backend):
        be.set_backend(backend)
        proto = actions.Prototype(frequency=2, start=1, stop=3)
        u = MagicMock()

        u.time = 0.5
        assert not proto.will_run(0, u)
        u.time = 1.5
        assert proto.will_run(2, u)
        u.time = 2.5
        assert not proto.will_run(1, u)
        u.time = 3.5
        assert not proto.will_run(4, u)

    def test_call(self, backend):
        be.set_backend(backend)
        proto = actions.Prototype()
        proto.will_run = MagicMock(return_value=True)
        proto._doit = MagicMock()
        u = MagicMock()

        proto(1, u)
        proto.will_run.assert_called_once_with(1, u)
        proto._doit.assert_called_once_with(1, u)


@pytest.mark.parametrize("backend", backends)
class TestBlowupCutoff:
    def test_above_cutoff(self, backend):
        be.set_backend(backend)
        blowup = actions.BlowupCutoff(cutoff=10)
        u = MagicMock()

        u.fields = [be.array([5, 9]), be.array([9, -6])]
        assert not blowup.above_cutoff(u)

        u.fields = [be.array([5, 10]), be.array([10, 20])]
        assert blowup.above_cutoff(u)

    def test_doit(self, backend):
        be.set_backend(backend)
        blowup = actions.BlowupCutoff(cutoff=10)
        u = MagicMock()

        u.fields = [be.array([5, 9]), be.array([9, -6])]
        blowup._doit(1, u)  # No exception should be raised

        u.fields = [be.array([5, 10]), be.array([10, 20])]
        with pytest.raises(Exception, match="Function values are above the cutoff"):
            blowup._doit(1, u)


@pytest.mark.parametrize("backend", backends)
class TestInfo:
    def test_doit(self, backend, caplog):
        be.set_backend(backend)
        info = actions.Info()
        u = MagicMock()
        u.time = 2.0

        with caplog.at_level(logging.INFO, logger="Info"):
            info._doit(3, u)
        assert "Iteration: 3, Time: 2.000000" in caplog.text
