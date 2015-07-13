# Discrete L2 Wasserstein Barycenters

Julia code (tested on OSX 10.10.4 and Julia v0.3) for constructing the discrete Wasserstein barycenters presented in the paper titled `Discrete Wasserstein Barycenters: Optimal Transport for Discrete Data` by Ethan Anderes, Steffen Borgwardt and Jacob Miller.

To excite the code in this repo you can first copy/clone all the source code by typing the followign command in the terminal.

```
$ git clone https://github.com/EthanAnderes/WassersteinBarycenterCode.
```

 Now simply move into the source directory and launch with the following commands.

```
$ cd WassersteinBarycenterCode/
$ julia
julia> Pkg.add("DataFrames")
julia> Pkg.add("PyCall")
julia> Pkg.add("PyPlot")
julia> Pkg.add("JuMP")
julia> Pkg.add("Clp")
julia> include("scripts/HeatDemand.jl")
```

Note, if you already have the following packges installed, then adding the packges with Pkg.add is not necssary. The packages `PyCall` and `PyPlot` allow Julia to call Python for generating figures. Please refer to the pacakge documentation for both `PyCall` and `PyPlot` for appropriate configuration.
