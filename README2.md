# UlamMethod.jl


<!-- [!["Documentation link"](src/assets/docu.svg)](https://70gage70.github.io/UlamMethod.jl/) -->

## Introduction

This package is an implementation of Ulam's method [^1] [^2] (see also Galerkin projection [^3]) for the discretization of a stochastic operator using pure Julia. Given a set of two-dimensional, one-step trajectories 
```math
(x_{0, 1}, y_{0, 1}) \to  (x_{T, 1}, y_{T, 1}), (x_{0, 2}, y_{0, 2}) \to  (x_{T, 2}, y_{T, 2}) \dots
```
defined in a domain, the essential goal of Ulam's method is to partition the domain into a series of non-intersecting regions and construct a transition probability matrix $P$ on these regions. In UlamMethod.jl, this is accomplished in two main steps

1. The user provides a geometry containing the data, and covering of the the domain is generated according to one of several different binning algorithms (lines, rectangles, and [Voronoi cells](https://en.wikipedia.org/wiki/Voronoi_diagram).)

2. The number of trajectories beginning in polygon $i$ and ending in polygon $j$ is used to create the entry $P_{i, j}$ of $P$ such that the edges of the domain are handled by a stochasticization algorithm which re-injects the data either according to interior-exterior trajectories or at a user-defined location.

The geometries which form the covering and the transition probability matrix are the main outputs.

## Installation

This package is in the Julia General Registry. In the Julia REPL, run the following code and follow the prompts:

```julia
using Pkg
Pkg.add("UlamMethod")
```

Access the functionality of the package in your code by including the following line:

```julia
using UlamMethod
```

## Quickstart

The core functionality is provided by 
```julia
ulam_method(traj, boundary, binner; reinj_algo)
``` 
where

- `traj`: A [`Trajectories`](@ref) object, holding the short-range trajectory data.
- `boundary`: A [`Boundary`](@ref) object, holding the geometry that defines the computational boundary.
- `binner`: A [`BinningAlgorithm`](@ref) that specifies the algorithm used to partition the boundary into bins.
- `reinj_algo`: A [`ReinjectionAlgorithm`](@ref) that specifies how trajectories pointing from nirvana to the interior should be reinjected. Default [`DataReinjection`](@ref).

Here are `10000` random trajectories in the domain $[0, 10]^2$

```julia
using UlamMethod

n_data = 10000
x0, xT = 10*rand(2, n_data), 10*rand(2, n_data)
traj = Trajectories(x0, xT)
```

We will take our domain to be the rectangular subset $[3, 5] \times [4, 8]$ and generate a covering with 40 rectangles.

```julia
xmin, xmax, ymin, ymax = 3, 5, 4, 8
boundary = Boundary(xmin, xmax, ymin, ymax)
binner = RectangleBinner(40)
```

Run Ulam's method.

```julia
ulam = ulam_method(traj, boundary, binner)
```

Use `P_closed(ulam)` to see the full transition probability matrix and `bins(ulam)` to see the list of bins.

## References

[^1]: Ulam, Stanislaw M. A collection of mathematical problems. No. 8. Interscience Publishers, 1960.

[^2]: Li, Tien-Yien. "Finite approximation for the Frobenius-Perron operator. A solution to Ulam's conjecture." Journal of Approximation theory 17.2 (1976): 177-186.

[^3]: Reddy, Junuthula Narasimha. Introduction to the finite element method. McGraw-Hill Education, 2019.

