# Advanced Usage

```@contents
Pages = ["advanced.md"]
Depth = 5
```

# High level overview

Use of this package proceeds as follows.

1. Create `Trajectories` using your trajectory data
2. Select the computational domain (`Boundary`). Data inside the boundary are considered "interior" and data outside the boundary are considered "exterior" or "in nirvana"[^1] [^2]. Use the `AutoBoundary` function to attempt to select a reasonable boundary automatically. In 2D, the `AutoBoundary2D` function is also available, which may find a tighter boundary than `AutoBoundary`.
3. Select a `BinningAlgorithm` and set any parameters it has.
4. Optionally, select a `ReinjectionAlgorithm` to handle the behavior of trajectories points from nirvana to the interior.
5. Run `ulam_method(traj, binner; reinj_algo)`.
6. Inspect and process the result with `P_closed`, `bins`, `membership` and other functions.
7. Visualize the result with `viz` and `viz!`.

# Loading trajectories

```@docs; canonical=false
Trajectories
```

# Defining Boundaries

```@docs; canonical=false
Boundary
```

# Binning Algorithms

```@docs; canonical=false
BinningAlgorithm
```

## 1D Algorithms

```@docs; canonical=false
LineBinner
```

## 2D Algorithms

```@docs; canonical=false
RectangleBinner
TriangleBinner
HexagonBinner
VoronoiBinner
```

## ≥3D Algorithms

```@docs; canonical=false
HyperRectangleBinner
```

# Reinjection algorithms

```@docs; canonical=false
ReinjectionAlgorithm
```

```@docs; canonical=false
DataReinjection
SourceReinjection
StationaryReinjection
```

# Computing the main results

```@docs; canonical=false
ulam_method
```

# Working with results

```@docs; canonical=false
UlamResult
P_closed
bins
bins_dis
membership
points
counts
```

# Visualization

```@docs; canonical=false
viz
viz!
```

# References

[^1]: Miron, Philippe, et al. "Transition paths of marine debris and the stability of the garbage patches." Chaos: An Interdisciplinary Journal of Nonlinear Science 31.3 (2021): 033101.

[^2]: In brief, nirvana is an extra state appended to an open system to close it; trajectories which point from inside the domain to the outisde of the domain transition to this nirvana state. Trajectories which point from outside the domain to the inside are transitions "from" nirvana - how exactly these data are reinjected is controlled by the `ReinjectionAlgorithm`.