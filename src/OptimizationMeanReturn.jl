using LinearAlgebra, ForwardDiff, Random, Plots

"""
    StochasticDominanceMeanReturn(ξ, ξ_0, p_ξ, p_ξ_0, SDorder; ε=1e-7, max_iter=50, n_particles=200, verbose=false, plot=false)

Solve a portfolio optimization problem under stochastic dominance constraints to maximize the expected return.

# Arguments
- `ξ::AbstractMatrix{<:Number}`: Portfolio return scenarios.
- `ξ_0::AbstractVector{<:Number}`: Benchmark return scenarios.
- `p_ξ::AbstractVector{<:Number}`: Probability distribution over portfolio scenarios.
- `p_ξ_0::AbstractVector{<:Number}`: Probability distribution over benchmark scenarios.
- `SDorder::Int`: Order of stochastic dominance.
- `ε::Float64` (optional, default = `1e-7`): Convergence tolerance.
- `max_iter::Int` (optional, default = `50`): Maximum iterations for optimization.
- `n_particles::Int` (optional, default = `200`): Number of particles for PSO.
- `verbose::Bool` (optional, default = `false`): Print convergence details.
- `plot::Bool` (optional, default = `false`): Plot the optimal asset allocation.

# Returns
- A tuple `(x_opt, t_opt)`, where:
  - `x_opt`: Optimized portfolio weights.
  - `t_opt`: Auxiliary optimization variable.

# Examples
```julia
StochasticDominanceMeanReturn([2 3; 5 6], [2,3,4], [0.5,0.5], [0.2,0.3,0.5], 2, plot=true)
```
# Solves the portfolio optimization problem under stochastic dominance constraints.
"""
function StochasticDominanceMeanReturn(ξ, ξ_0,p_ξ, p_ξ_0,SDorder;ε=1e-4, max_iter=50,n_particles = 200,verbose=false,plot=false)
    length_x = size(ξ)[1]
    
    # This is a psuedo initialization (For PSO)
    x= rand(length_x)
    t=rand()
    if SDorder == 1
        error("First order is not available")
        println("Try order: $SDorder+1")
    end
    p=SDorder-1.0
    # Define the combined function without passing tilde_t_1 explicitly
    function combined_fun(vars)
        x, t = unpack_vars(vars)
    
        # Compute constraint functions
        simplex_constraints = vcat(
            simplex1(x),
            simplex2(x)          
        )
        stochastic_dominance_constraints = vcat( g_p(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0),
        g_p_minus_1(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0),
        g_bar(x, ξ, ξ_0, p, p_ξ, p_ξ_0),
        boundary_moments(x, ξ, p_ξ,p_ξ_0,ξ_0, p))
        
        # Compute Mean Return function (objective)
        meanreturn = -MeanReturn(x,ξ,p_ξ)
      
       # Apply penalty
       constraint_penalty1 =  norm(100 .*simplex_constraints,Inf) 
       constraint_penalty2 =  norm(100 .*stochastic_dominance_constraints,Inf )
       
       if constraint_penalty1 < 1e-2
        constraint_penalty1=0.0
       end
       if constraint_penalty2 < 1e-2
        constraint_penalty2=0.0
       end
       result= meanreturn .+ constraint_penalty1 .+ constraint_penalty2 
       if result >0 
       result = 1000*result
       end      
       return result
    end
    function unpack_vars(vars)
        x = vars[1:length_x]                # Extract the first d elements for x
        t = vars[end]                # Last element is t
        return x,t
    end
    vars0 = vcat(x, t)
    result= PSO(combined_fun,vars0;n_particles = n_particles,maxEval = max_iter)
    
    x_opt, t_opt = unpack_vars(result.xMin)
    vars0 = vcat(x_opt,t_opt)
    
    ## Check Convergence
    
    ### If converged return else try newton 
    if CheckConvergenceSimplex(x_opt)+ CheckConvergenceSD(x_opt, t_opt, ξ, ξ_0, p, p_ξ, p_ξ_0) ≤ ε
        if verbose
            if CheckConvergenceSimplex(x_opt)+ CheckConvergenceSD(x_opt, t_opt, ξ, ξ_0, p, p_ξ, p_ξ_0) ≤ ε
                println("Converged")
                println("Simplex Constraints residuals: ", CheckConvergenceSimplex(x_opt))
                println("Stochastic Dominance Constraints residuals: ", CheckConvergenceSD(x_opt, t_opt, ξ, ξ_0, p, p_ξ, p_ξ_0))    
                else
                println("Not Converged:Consider increasing max_iter and n_particles")
                println("Simplex Constraints residuals: ", CheckConvergenceSimplex(x_opt))
                println("Stochastic Dominance Constraints residuals: ", CheckConvergenceSD(x_opt, t_opt, ξ, ξ_0, p, p_ξ, p_ξ_0))    
            end
        end
        if plot
            plotOptimalAssetAllocationMeanReturn(x_opt, ξ, p_ξ, ξ_0, p_ξ_0)
        end
      return x_opt, t_opt
    else
     function combined_gradient(vars)
           return ForwardDiff.gradient(combined_fun, vars)
     end
    # Main optimization loop
    vars_opt = vars0
       
    result = Newton(combined_fun, combined_gradient, vars_opt; maxEval=max_iter)
    x_opt,t_opt = unpack_vars(result.xMin)
        if verbose
            if CheckConvergenceSimplex(x_opt)+ CheckConvergenceSD(x_opt, t_opt, ξ, ξ_0, p, p_ξ, p_ξ_0) ≤ ε
                println("Converged")
                println("Simplex Constraints residuals: ", CheckConvergenceSimplex(x_opt))
                println("Stochastic Dominance Constraints residuals: ", CheckConvergenceSD(x_opt, t_opt, ξ, ξ_0, p, p_ξ, p_ξ_0))    
                else
                println("Not Converged:Consider increasing max_iter and n_particles")
                println("Simplex Constraints residuals: ", CheckConvergenceSimplex(x_opt))
                println("Stochastic Dominance Constraints residuals: ", CheckConvergenceSD(x_opt, t_opt, ξ, ξ_0, p, p_ξ, p_ξ_0))    
            end
        end
        if plot
            plotOptimalAssetAllocationMeanReturn(x_opt, ξ, p_ξ, ξ_0, p_ξ_0)
        end
      return  x_opt, t_opt
    end
end

