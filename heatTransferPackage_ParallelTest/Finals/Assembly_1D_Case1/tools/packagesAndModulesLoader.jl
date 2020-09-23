include("../modules/FDMD1HeatTransfer.jl")
include("../modules/LBM_heat_D1Q2_v3.jl")
@info("Loading modules")
using .FDMHeatTransfer
using .D1Q2
using Pkg
@info("Loading modules: successfull")

if haskey(Pkg.installed(),"Plots")
    @info("Loading Plots.jl: ")
    using Plots
    @info("Loading Plots.jl: successfull")
else
    @error("You have not installed Plots.jl Package.")
end

if haskey(Pkg.installed(),"BenchmarkTools")
    @info("Loading BenchmarkTools.jl: ")
    using BenchmarkTools
    @info("Loading BenchmarkTools.jl: successfull")
else 
    @error("You have not installed BenchmarkTools.jl Package.")
end
# if haskey(Pkg.installed(),"ArgParse")
#     @info("Loading ArgParse.jl: ")
#     using ArgParse
#     @info("Loading ArgParse.jl: successfull")
# else 
#     @error("You have not installed ArgParse.jl Package.")
# end

@info("Compiling ...")
