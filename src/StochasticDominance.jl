module StochasticDominance

# Import necessary Julia libraries
using LinearAlgebra, ForwardDiff, Random, Plots, GR


# Include all submodules
include("AuxiliaryFunc.jl")
include("NewtonOpt.jl")
include("PSOOpt.jl")
include("OptimizationMeanReturn.jl")
include("OptimizationRiskFunction.jl")

# Import submodules
using .AuxiliaryFunc, .NewtonOpt, .PSOOpt

# Export functions from main functions
export StochasticDominanceMeanReturn, StochasticDominanceRiskMeasure


# Export functions from AuxiliaryFunc
export safe_exponent, simplex1, simplex2, g_p_minus_1, g_p, g_bar, boundary_moments, 
       CheckConvergenceSimplex, CheckConvergenceSD, MeanReturn, 
       plotOptimalAssetAllocationMeanReturn, RiskFunction, BenchmarkRiskFunction, 
       plotOptimalAssetAllocationRiskFunction

# Export optimization solvers
export Newton, PSO

end
