# Theory and Implementation

UNDER CONSTRUCTION.

## Calculation Hierarchy
- The user provides `UlamTrajectories` via file or manual input. This type contains `x0`, `y0`, `xT`, `yT`.
- The user provides an `UlamDomain` which contains the `domain`, an `UlamPolygon` defining the boundary of nirvana. Also provided is the type and number of polygons requested, and the type and location of the reinjection algorithm.
- The highest level function is `ulam_method` which calls the various steps in the calculation.
- First, `ulam_binner` is used. This selects and runs the appropriate binning algorithm based on the user's choice. Every binning algorithm acts on the bounding box of `domain`. The output of every binning algorithm is vector of `UlamPolygon`s which represent a covering of the bounding box.
- Afterwards, the `domain` is intersected with the binning result, again giving a vector of `UlamPolygon` but such that every data point is now inside `domain`. The intersection is such that if the result of an intersection is a multipolygon, the polygon with the largest area is taken.
- The result of `ulam_binner` "actual" Ulam's method calculation, which is contained in `ulam_nirvana`. This function takes the trajectories and polygons and constructs the transition matrix defining the Markov chain on states representing the polygons. Polygons which contain no `(x0, y0)` points are discarded. The largest strongly connected component of the transition matrix is extracted. Both of these "cleaning" operations have the effect of deleting trajectories which start or end in any of the removed polygons.