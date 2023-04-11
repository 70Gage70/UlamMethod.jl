# Theory and Implementation

UNDER CONSTRUCTION.

## Calculation Hierarchy
- The user provides UlamTracjectories via file or manual input. This type contains `x0`, `y0`, `xT`, `yT`.
- The user provides an UlamDomain. The minimum specification of such a domain is the corners. Also optionally included
is a `domain` which is a further restriction that the user can provide, namely that all the data should lie
inside the polygon defined by the domain. Also provided is the type and number of polygons requested, and the
type and location of the reinjection algorithm.
- The highest level function is `ulam_method` which calls the various steps in the calculation.
- First, `ulam_binner` is used. This selects and runs the appropriate binning algorithm based on the user's choice. The output of every binning algorithm is vector of `UlamPolygon`s which represent a covering of the box defined by the corners.
- Afterwards, if the user provided a `domain`, this domain is intersected with the binning result, again giving a vector of `UlamPolygon` but such that every data point is now inside `domain`. Note that the result of this is that any points inside the corners but outside the domain are treated as belonging to nirvana.
- The result of this is then passed to the "actual" Ulam's method calculation, which is contained in `ulam_method`. This function takes the trajectories and polygons and constructs the transition matrix defining the Markov chain on states representing the polygons. Polygons which contain no `(x0, y0)` points are discarded. The largest strongly connected component of the transition matrix is extracted. Both of these "cleaning" operations have the effect of deleting trajectories which start or end in any of the removed polygons.