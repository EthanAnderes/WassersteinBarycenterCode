# Discrete L2 Wasserstein Barycenters

This repo contains Julia code for constructing the discrete Wasserstein barycenters presented in the paper titled `Discrete Wasserstein Barycenters: Optimal Transport for Discrete Data` by Ethan Anderes, Steffen Borgwardt and Jacob Miller (a local copy of the paper can be found in directory `paper`).


## Version information

All the code was tested on Julia 0.4 (installation instructions can be found at [julialang.org](http://julialang.org)). The exact version information for the system which generated the figures are given as follows. 
```julia
julia> versioninfo()
Julia Version 0.4.6-pre+4
Commit 03c072d (2016-03-19 15:00 UTC)
Platform Info:
  System: Darwin (x86_64-apple-darwin15.4.0)
  CPU: Intel(R) Core(TM) i7-4850HQ CPU @ 2.30GHz
  WORD_SIZE: 64
  BLAS: libopenblas (USE64BITINT DYNAMIC_ARCH NO_AFFINITY Haswell)
  LAPACK: libopenblas64_
  LIBM: libopenlibm
  LLVM: libLLVM-3.3

julia> Pkg.installed("PyPlot")
v"2.1.1"

julia> Pkg.installed("JuMP")
v"0.12.2"

julia> Pkg.installed("Clp")
v"0.2.0"

julia> Pkg.installed("PyCall")
v"1.4.0"

julia> Pkg.installed("PyPlot")
v"2.1.1"

```
This information is intended to enable anyone to recreate the exact system which generated the figures in the paper. That said, we expect the code to be forward compatible with stable releases.


*Note:* Most of the code to compute the Wasserstein barycenters is done in Julia. However, the actual plotting of the figures figures is done using python modules modules `matplotlib` and `mpl_toolkits.basemap` (calling from within Julia via the packages `PyPlot` and `PyCall`). Therefore, if you want to generate the graphics, you will need to have these modules working in python.




## Running the code to generate the figures

To execute the code in this repo you can first copy/clone all the source code by typing the following command in the terminal.

```bash
$ git clone https://github.com/EthanAnderes/WassersteinBarycenterCode
```

Now move to the repo directory and launch Julia.

```bash
$ cd WassersteinBarycenterCode/
$ julia
```

Now add the required packages and launch the script file `HeatDemand.jl`

```julia
julia> Pkg.add("DataFrames")
julia> Pkg.add("PyCall")
julia> Pkg.add("PyPlot")
julia> Pkg.add("JuMP")
julia> Pkg.add("Clp")
julia> include("scripts/HeatDemand.jl")
```

Note, if you already have the following packages installed, then adding the packages with Pkg.add is not necessary. The packages `PyCall` and `PyPlot` allow Julia to call Python for generating figures. Please refer to the package documentation for both `PyCall` and `PyPlot` for appropriate configuration.
