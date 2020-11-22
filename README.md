# SnaggedNecklace

Constrained path integral Monte Carlo (PIMC) and numerical matrix multiplication (NMM) for the computation of the quantum potential of mean force (PMF) of tiny systems.

Tested with Julia 1.5.


## Installation

```
pkg> add https://github.com/0/SnaggedNecklace.jl.git
```

In order to run the example scripts in `examples/`, you will also need to
```
pkg> add ArgParse
```

### Application project

If you're working with a clone of this repository, you can use the basic application project in `examples/`, which already has both `SnaggedNecklace` and `ArgParse` as dependencies.
From the repository root, run
```
julia --project=examples
```
and then
```
pkg> dev .
```
to create `examples/Manifest.toml` with a development version of `SnaggedNecklace`.


## References

The binning method is from: Ambegaokar, V., & Troyer, M. (2010). Estimating errors reliably in Monte Carlo simulations of the Ehrenfest model. American Journal of Physics, 78(2), 150-157. [doi:10.1119/1.3247985](https://aapt.scitation.org/doi/abs/10.1119/1.3247985), [arXiv:0906.0943](https://arxiv.org/abs/0906.0943).


## Acknowledgements

Thanks to Kevin Bishop for helping to verify this implementation!


## License

Provided under the terms of the MIT license.
See `LICENSE` for more information.
