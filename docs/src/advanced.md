# Advanced Usage

## [Binning Algorithms](@id binning)

## [Stochasticization Algorithms](@id stoc)

## Under construction

## Extra

- Bin the data according to the covering to generate a transition probability matrix such that any data outside the computational domain is lumped into a so-called "nirvana" state [^4]. Trajectories leaving the domain are reinjected using one of two possible [algorithms](#reinjection-algorithms).
- Transition path theory statistics [^5] are computed on the stationary, time-homogenous Markov chain induced by the transition probability matrix.
- The outputs are written to a [.h5](https://github.com/JuliaIO/HDF5.jl) file.

Prepare the data into one of two formats: [.mat](https://github.com/JuliaIO/MAT.jl) or [.h5](https://github.com/JuliaIO/HDF5.jl) such that the head of the file contains four variables, `x0`, `xT`, `y0` and `yT`. See `test/x0x5-NA-undrogued.mat` and `test/x0x5-NA-undrogued.h5` for example trajectory data from undrogued drifters in the North Atlantic obtained from the NOAA GDP [^6].