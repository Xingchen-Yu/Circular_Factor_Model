"""Circular (spherical) latent factor model for roll-call data, in JAX.

Python/JAX port of Yu & Rodriguez (2021), "Spatial voting models in circular
spaces: A case study of the U.S. House of Representatives", AOAS 15(4).
"""

import os

# WSL/consumer GPUs often cannot satisfy XLA's default 75% preallocation.
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
