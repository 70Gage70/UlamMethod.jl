# UlamMethod.jl

```@contents
Pages = ["index.md"]
Depth = 5
```

# Introduction

This is the documentation for the Julia package [`UlamMethod.jl`](https://github.com/70Gage70/UlamMethod.jl).

This package is an implementation of Ulam's method [^1] [^2] (see also Galerkin projection [^3]) for the discretization of a stochastic operator using pure Julia. Given a set of one-step trajectories 
```math
\mathbf{x}_{0, 1} \to  \mathbf{x}_{T, 1}, \mathbf{x}_{0, 2} \to  \mathbf{x}_{T, 2} \dots
```
defined in a subset of $\mathbb{R}^N$, the essential goal of Ulam's method is to partition the domain into a series of non-intersecting regions and construct a transition probability matrix $P$ on these regions.  

## Features
- Built on the [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl) computational geometry framework.
- Supports trajectory data in arbitrary dimensions.
- Supports automatic boundary construction.
- Multiple 2D algorithms for partitioning to triangles, rectangles, hexagons and adaptively sized [Voronoi cells](https://en.wikipedia.org/wiki/Voronoi_diagram).
- Multiple stochasticization algorithms.

# Installation

This package is in the Julia General Registry. In the Julia REPL, run the following code and follow the prompts:

```julia
import Pkg
Pkg.add("UlamMethod")
```

Access the functionality of the package in your code by including the following line:

```julia
using UlamMethod
```

## Quickstart

The core functionality is provided by 
```julia
ulam_method(traj, binner; reinj_algo)
``` 
where

- `traj`: A `Trajectories` object, holding the short-range trajectory data.
- `boundary`: A `Boundary` object, holding the geometry that defines the computational boundary.
- `reinj_algo`: A `ReinjectionAlgorithm` that specifies how trajectories pointing from nirvana[^4] to the interior should be reinjected. Default [`DataReinjection`](@ref).

Here are `10000` random trajectories in the domain $[0, 10]^2$

```julia
using UlamMethod
import Random; Random.seed!(1234) # reproducible randomness

n_data = 10000
x0, xT = 10*rand(2, n_data), 10*rand(2, n_data)
traj = Trajectories(x0, xT)
```

We will take our domain to be the rectangular subset $[3, 5] \times [4, 8]$ and generate a covering with 40 rectangles. This covering is defined inside a `Boundary` object, which can be quickly created in 2D using the syntax `Boundary(xmin, xmax, ymin, ymax)`. We then call `ulam_method` to run the main calculation.

```julia
xmin, xmax, ymin, ymax = 3, 5, 4, 8
boundary = Boundary(xmin, xmax, ymin, ymax)
binner = RectangleBinner(40, boundary)

ulam = ulam_method(traj, binner)
```

- `P_closed(ulam)` gives the full transition probability matrix.
-  `bins(ulam)` gives the bins and `points(bins(ulam))` gives their vertices.
- `membership(points, ulam)` returns the bin membership of a `Dim x n_points` matrix `points`.

# Citation

!!! note

    Please use the following citation if you use this package in your research.

    ```
    @article{bonner2023improving,
    title={Improving the stability of temporal statistics in transition path theory with sparse data},
    author={Bonner, Gage and Beron-Vera, FJ and Olascoaga, MJ},
    journal={Chaos: An Interdisciplinary Journal of Nonlinear Science},
    volume={33},
    number={6},
    year={2023},
    publisher={AIP Publishing}
    }
    ```

Initial development of this package was supported by the National Science Foundation.

# References

[^1]: Ulam, Stanislaw M. A collection of mathematical problems. No. 8. Interscience Publishers, 1960.

[^2]: Li, Tien-Yien. "Finite approximation for the Frobenius-Perron operator. A solution to Ulam's conjecture." Journal of Approximation theory 17.2 (1976): 177-186.

[^3]: Reddy, Junuthula Narasimha. Introduction to the finite element method. McGraw-Hill Education, 2019.

[^4]: In brief, nirvana is an extra state appended to an open system to close it; trajectories which point from inside the domain to the outisde of the domain transition to this nirvana state. Trajectories which point from outside the domain to the inside are transitions "from" nirvana - how exactly these data are reinjected is controlled by the `ReinjectionAlgorithm`.
