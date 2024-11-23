# Compiling UlamMethod.jl Into a Standalone App

## Creating The App

This experimental feature allows you to compile `UlamMethod.jl` into a standalone executable that can be used independently of Julia. Currently, only 2D data are supported. Note that creating the app requires a julia installation.

1. Create a folder where you want the app to be stored and `cd` into it. For this example, we'll use `UlamApp`.

2. Clone [https://github.com/70Gage70/UlamMethod.jl](https://github.com/70Gage70/UlamMethod.jl) into `UlamApp`. If you use [GitHub CLI](https://cli.github.com/) this is accomplished quickly by: `gh repo clone 70Gage70/UlamMethod.jl`. 

3. Run the following command from the terminal. This will compile the app, taking 5-10 minutes. It will also add the package [`PackageCompiler`](https://github.com/JuliaLang/PackageCompiler.jl/tree/master) to your main environment. Run `julia -e 'import Pkg; Pkg.rm("PackageCompiler")` afterwards if you don't want it any more.
```
julia -e 'import Pkg; Pkg.add("PackageCompiler"); Pkg.develop(path = "UlamMethod.jl"); using PackageCompiler, UlamMethod; create_app("UlamMethod.jl/", "UlamMethodCompiled")'
```
The path to the executable is now `.../UlamApp/UlamMethodCompiled/bin/UlamMethod`. This can e.g. be added to your PATH for convenient access.

## Using The App

The expected syntax is `./UlamMethod ARG1 ARG2 ... ARG9`, where each `ARG` is defined as follows

- `ARG1`: The path to the file containing the input data. This expects a [`.mat`](https://github.com/JuliaIO/MAT.jl) file that contains the variables `"x0", "xT", "y0"` and `"yT"`.
- `ARG2`: The file to write the output. This should also be a `.mat` file.
- `ARG3`: The type of binner. Options are `[rec, tri, hex, vor]`.
- `ARG4`: The number of bins.
- `ARG5 - ARG8`: These args control the fraction of data taken to be in nirvana on each side of the rectangular boundary in the order `left, right, bottom, top`. 
- `ARG9`: The reinjection algorithm. Options are `[data, stat, unif]`.


### Example

The following terminal command will generate some test data in the file `traj-test.mat`:

```
julia -e 'import Pkg; Pkg.activate(temp = true); Pkg.add(["UlamMethod", "MAT"]); using UlamMethod, MAT; traj = Trajectories(2, 1000); matwrite("traj-test.mat", Dict("x0" => traj.x0[1,:], "y0" => traj.x0[2,:], "xT" => traj.xT[1,:],"yT" => traj.xT[2,:]))'
```

We'll apply Ulam's method with 100 rectangular bins and `"data"` reinjection (the default). We'll let `2%` of the data be in nirvana on the left side of the boundary. The output file will be `ulam-out.mat`. Due to precompilation issues, the very first time you run the app, it will be slow.

```
./UlamMethod traj-test.mat ulam-out.mat rec 100 0.02 0.0 0.0 0.0 data
```

The output `.mat` file contains the transition probability matrix `"P"` as well as the vertices defining the associated bins `"bins_verts"`. Refer to the `"info"` field for more information regarding the format of the output.