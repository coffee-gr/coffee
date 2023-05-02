import pytest

from coffee.backend import Backend


class MockBackend:
    def __init__(self):
        self.calls = []

    def __getattr__(self, name):
        def method(*args, **kwargs):
            self.calls.append((name, args, kwargs))
            return f"{name} called"

        return method


class TestBackend:
    @pytest.fixture
    def backend(self):
        return Backend()

    def test_set_backend_default(self, backend):
        backend.set_backend()
        assert backend.backend_name == "numpy"

    def test_set_backend_unsupported(self, backend):
        with pytest.raises(ValueError):
            backend.set_backend("unsupported")

    def test_get_backend_unsupported(self, backend):
        with pytest.raises(ValueError):
            backend._get_backend("unsupported")

    def test_backend_methods(self, backend, monkeypatch):
        mock_backend = MockBackend()
        monkeypatch.setattr(backend, "_backend", mock_backend)

        backend.abs(1)
        backend.array([1, 2, 3])
        backend.linspace(0, 1, 5)
        backend.zeros((3, 3))
        backend.sum([1, 2, 3])

        expected_calls = [
            ("abs", (1,), {}),
            ("array", ([1, 2, 3],), {}),
            ("linspace", (0, 1, 5), {}),
            ("zeros", ((3, 3),), {}),
            ("sum", ([1, 2, 3],), {}),
        ]

        assert mock_backend.calls == expected_calls
