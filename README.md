# Circular Factor Model

[![DOI](https://zenodo.org/badge/271097395.svg)](https://zenodo.org/badge/latestdoi/271097395)

Code for

> Yu, X. & Rodríguez, A. (2021). *Spatial voting models in circular spaces: A
> case study of the U.S. House of Representatives.* Annals of Applied
> Statistics 15(4), 1897–1922.

| Folder | Contents |
|---|---|
| [`R/`](R/) | Original R + Rcpp implementation used for the paper (no longer maintained), with the 100th–116th House `.ord` data. See [`R/README.md`](R/README.md). |
| [`python/`](python/) | **Default.** JAX port (`circfactor`) that runs on a GPU, tested against high-precision references and validated against the original 116th House results. See [`python/README.md`](python/README.md). |

The model is also available as an R package:
[SLFM1D](https://github.com/Xingchen-Yu/SLFM1D).
