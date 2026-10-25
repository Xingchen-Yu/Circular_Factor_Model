import circfactor  # noqa: F401  (sets XLA env defaults before jax import)
import jax

# Tests run in float64 by default; float32 paths are exercised with explicit dtypes.
jax.config.update("jax_enable_x64", True)
