# Circular Factor Model

[![DOI](https://zenodo.org/badge/271097395.svg)](https://zenodo.org/badge/latestdoi/271097395)

Code for

> Yu, X. & Rodríguez, A. (2021). [*Spatial voting models in circular spaces: A
> case study of the U.S. House of Representatives.*](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-15/issue-4/Spatial-voting-models-in-circular-spaces--A-case-study/10.1214/21-AOAS1454.full)
> Annals of Applied Statistics 15(4), 1897–1922.

Explore the results for every House since 1987 at **[fullcirclecongress.org](https://fullcirclecongress.org)**
(source in [`website/`](website/)).

| Folder | Contents |
|---|---|
| [`python/`](python/) | JAX implementation (`circfactor`) that runs on a GPU, tested against high-precision references and validated against the original 116th House results. See [`python/README.md`](python/README.md). |

## Deprecated

These are kept for reference only and are no longer maintained:

* [`R/`](R/): the original R + Rcpp implementation used for the paper, with the 100th–116th House
  `.ord` data. See [`R/README.md`](R/README.md).
* [SLFM1D](https://github.com/Xingchen-Yu/SLFM1D): an earlier R package of the model.

Use the Python implementation above for new work.
