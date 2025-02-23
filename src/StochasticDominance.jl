module StochasticDominance

# Import necessary Julia libraries
using LinearAlgebra, ForwardDiff, Random, Plots, GR


# Include all submodules
include("AuxiliaryFunc.jl")
include("NewtonOpt.jl")
include("PSOOpt.jl")
include("OptimizationMeanReturn.jl")
include("OptimizationRiskFunction.jl")
include("OptimizationVerifyDominance.jl")

# Import submodules
using .AuxiliaryFunc, .NewtonOpt, .PSOOpt

# Export functions from main functions
export optimize_max_return_SD, optimize_min_riskreturn_SD, optimize_min_risk,verify_dominance


# Export functions from AuxiliaryFunc
export safe_exponent, simplex1, simplex2, g_p_minus_1, g_p, g_bar, boundary_moments, 
       convergence_simplex, convergence_SD, expected_portfolio_return, 
       plot_max_return, riskfunction_asset_allocation, riskfunction, 
       plot_min_return,randomVariable

# Export optimization solvers
export Newton, PSO

end
