__all__ = ['diffusion', 'poisson', 'schrodinger', 'schrodinger_poisson']

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__)))

from schrodinger import solve_schrodinger
from poisson import solve_poisson
from diffusion import solve_diffusion
from schrodinger_poisson import *
