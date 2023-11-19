# Import standard code base and initialize NumPy backend
from coffee.settings import init
from coffee.backend import backend as be
my_backend = be.set_backend("jax")
init(my_backend)

from coffee.diffop import fd

diffop = fd.FD12()

# Set up grid and f(x) = sin(x)
N  = 101
x  = be.linspace(-1,1,N)
dx = x[1] - x[0]
fn = be.sin(x)

# Take a derivative
dfn = diffop(fn, dx)

print(dfn - be.cos(x))