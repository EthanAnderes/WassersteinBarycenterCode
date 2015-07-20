# Discrete L2 Wasserstein Barycenters

Julia code (tested on OSX 10.10.4 and Julia v0.3) for constructing the discrete Wasserstein barycenters presented in the paper titled `Discrete Wasserstein Barycenters: Optimal Transport for Discrete Data` by Ethan Anderes, Steffen Borgwardt and Jacob Miller (a local copy of the paper can be found in directory `paper`).

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
