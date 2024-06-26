# Meshes Formulation

- A `Boundary` is a `Meshes.Mesh`. For example, a rectangular region:


```julia
points = [(0,0),(1,0),(0,1),(1,1)]
connec = connect.([(1,2,4,3)], Ngon)
mesh = SimpleMesh(points, connec)
```
- The user should provide a `Boundary`, then all data outside the boundary are considered to be in nirvana. 

- The user should provide a `AbstractBinner` which calculates a `SimpleMesh` from a `Boundary`.

