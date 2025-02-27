using LinearAlgebra, ForwardDiff, Random, Plots,Distributions

"""
    optimize_min_riskreturn_SD(ξ, ξ_0,SDorder; β=0.5,r=2.0,p_ξ, p_ξ_0, ε=1e-5, plot=false, verbose=false, max_ipot=50, max_iter=50, n_particles=50)

    Solve a portfolio optimization problem under stochastic dominance risk constraints

# Arguments
- `ξ::AbstractMatrix{<:Number}`: Portfolio return scenarios.
- `ξ_0::AbstractVector{<:Number}`: Benchmark return scenarios.
- `SDorder::S`: Order of stochastic dominance.
- `r::Number` (optional, default = `2.0`): norm of risk measure.
- `β::Number` (optional, default = `0.5`): Confidence level for risk measure.
- `p_ξ::AbstractVector{<:Number}`: Probability distribution over portfolio scenarios.
- `p_ξ_0::AbstractVector{<:Number}`: Probability distribution over benchmark scenarios.
- `ε::Float64` (optional, default = `1e-7`): Convergence tolerance.
- `plot::Bool` (optional, default = `false`): Plot the optimal asset allocation.
- `verbose::Bool` (optional, default = `false`): Print convergence details.
- `max_ipot::Integer` (optional, default = `50`): Maximum number of attempts to push the objective towards maximization. (Interior point optimization)
- `max_iter::Integer` (optional, default = `50`): Maximum iterations for optimization.
- `n_particles::Integer` (optional, default = `200`): Number of particles for PSO.



# Returns
- A tuple `(x_opt, q_opt, t_opt)`, where:
  - `x_opt`: Optimized portfolio weights.
  - `q_opt`: Optimized risk threshold parameter.
  - `t_opt`: Auxiliary optimization variable.

# Examples
```julia
# Define portfolio return scenarios (rows: assets, columns: scenarios)
ξ = [0.02  0.05 -0.01;  # Asset 1 returns under different scenarios
     0.03  0.06  0.02]   # Asset 2 returns under different scenarios

# Define benchmark return scenarios (column vector)
ξ_0 = [0.01, 0.04, 0.00]  # Benchmark returns in the same scenarios

# Define probability distributions for portfolio and benchmark scenarios
p_ξ = [0.3, 0.4, 0.3]     # Probabilities for each scenario
p_ξ_0 = [0.4, 0.4, 0.2]   # Probabilities for benchmark scenarios

# Set the order of stochastic dominance
SDorder = 2  # Second-order stochastic dominance

# Set the risk parameter
β=0.5

# Set the norm of the risk measure
r=2.0

x_opt, q_opt,t_opt= optimize_min_riskreturn_SD(ξ, ξ_0,SDorder;p_ξ, p_ξ_0,β,r,ε=1e-6,verbose=true,plot=true)
```
# Solves the portfolio optimization problem under stochastic dominance constraints.
"""
function optimize_min_riskreturn_SD(ξ::AbstractMatrix{T}, 
    ξ_0::AbstractVector{T}, 
    SDorder::S; 
    β::T=0.5, 
    r::T=2.0, 
    p_ξ=nothing, 
    p_ξ_0=nothing, 
    plot::Bool=false, 
    verbose::Bool=false, 
    ε::Float64=1e-5, 
    max_ipot::Integer=20, 
    max_iter::Integer=50, 
    n_particles::Integer=50) where {T<:Number, S<:Number}
        
    # Validate inputs
    if SDorder == 1
        error("First order is not available. Try order: $(SDorder+1)")
    elseif size(ξ, 2) ≠ length(ξ_0)
        error("The length of portfolio and benchmark scenarios are not equal. They must be equal to make the stochastic dominance comparison.")
    end
    # This is a psuedo initialization (For PSO)
    max_ipot = max_ipot  # Maximum times IP Optimization tries
    
    length_x = size(ξ, 1) # Lenth of portfolio
    max_ipot = max_ipot  # Maximum times IP Optimization tries
    p = SDorder - 1.0  # Define the dominance order

    # PSO Initialization
    x= rand(length_x)
    t=rand()    

    # Handling probabilities of ξ_0 and x'*ξ
    if p_ξ_0 !== nothing
        if length(p_ξ_0) == length(ξ_0)
            X = randomVariable(ξ_0, p_ξ_0)
        else
            error("The length of benchmark scenarios and their respective probabilities are not equal. They must be equal.")
        end
    else
        X = randomVariable(ξ_0)
    end
    p_ξ_0= X.probabilities
    if p_ξ !== nothing
        if length(p_ξ) == length(vec(x' * ξ))
            Y = randomVariable(vec(x' * ξ), p_ξ)
        else
            error("The length of portfolio scenarios and their respective probabilities are not equal. They must be equal.")
        end
    else
        Y = randomVariable(vec(x' * ξ))
    end
    p_ξ= Y.probabilities
    
    # Function for optimization (Simplex constraints (two) and stochastic dominance constraints of given order (Four))
    function combined_fun(vars)
        x, t = unpack_vars(vars)
           constraints = vcat(
                simplex1(x),
                simplex2(x),
                g_p(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0),
                g_p_minus_1(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0),
                g_bar(x, ξ, ξ_0, p, p_ξ, p_ξ_0),
                boundary_moments(x, ξ, p_ξ, p_ξ_0, ξ_0, p))
            result = norm(constraints)
        return result
    end
    
    # Function to unpack variables from optimization vector
    function unpack_vars(vars)
        x = vars[1:length_x]                
        t = vars[end]                
        return x, t
    end

    # Gradient of the combined function (for Newton's method)
    function combined_gradient(vars)
        return ForwardDiff.gradient(combined_fun, vars)
    end

    # Initialize optimization variables
    vars0 = vcat(x, t)
    # Optimize using PSO (initial) to find atleast one feasible point
    result = PSO(combined_fun, vars0; n_particles= n_particles, maxEval=max_iter,ε=ε)
    # Result after first try 
    x_opt, t_opt = unpack_vars(result.xMin)
    x_opt_best = x_opt; t_opt_best = t_opt
    q_opt_best = randn()

    # Trying 10 times using PSO to find at least one feasible point
    failure_count = 0 
    while true
        vars0 = vcat(x_opt_best, t_opt_best);
        result = PSO(combined_fun, vars0; n_particles=n_particles, maxEval=max_iter,Initialization=true,ε=ε)
        x_opt, t_opt = unpack_vars(result.xMin)
    
        if convergence_simplex(x_opt_best) + convergence_SD(x_opt_best, t_opt_best, ξ, ξ_0, p, p_ξ, p_ξ_0) <  ε
                break # Break the outer loop if a feasible point is already found
        elseif  convergence_simplex(x_opt) + convergence_SD(x_opt, t_opt, ξ, ξ_0, p, p_ξ, p_ξ_0) < convergence_simplex(x_opt_best) + convergence_SD(x_opt_best, t_opt_best, ξ, ξ_0, p, p_ξ, p_ξ_0) 
            x_opt_best = x_opt
            t_opt_best = t_opt # Update the best optimal points (x and t) if a better one is found
            failure_count = 0
        else
            failure_count += 1
            if failure_count >= max_ipot
                break  # Break the outer while loop after max_ipot consecutive failures
            end
        end 
    end 

    # If x^Tξ dominates ξ_0 in the given SDorder (ConvergenceRate_one_feasible_point < ε), otherwise dominance is not possible.
    ConvergenceRate_one_feasible_point = convergence_simplex(x_opt_best) + convergence_SD(x_opt_best, t_opt_best, ξ, ξ_0, p, p_ξ, p_ξ_0)

    # Find the optimal q (Risk function) only if at least one feasible point is found
    if ConvergenceRate_one_feasible_point < ε
    q_opt = optimize_min_risk(vec(x_opt_best'*ξ),p_ξ;β=β,r=r,ε=1e-5,max_iter=50,n_particles=1000,max_ipot=50) 
    q_opt_best = q_opt
    end
    failure_count = 0    # Reset failure counter
        #	╭────────────────────────────────────────────────────────────────
        #   │	Interior point optimization
    while true
            q_opt_best
            vars0 = vcat(x_opt_best,t_opt_best);
            result = PSO(combined_fun, vars0; n_particles= Int((n_particles)*0.5), maxEval=max_iter,Initialization=true,ε=ε)
            x_opt, t_opt = unpack_vars(result.xMin)
            
            # If solution is still infeasible, perturb variables and reattempt using Newton’s method
            if (convergence_simplex(x_opt) + convergence_SD(x_opt, t_opt, ξ, ξ_0, p, p_ξ, p_ξ_0) > ε) && p>1
                a=rand()
                temp_x = x_opt_best.+ rand(Uniform(-a,a),length_x) 
                vars0 = vcat(abs.(temp_x)/sum(abs.(temp_x)),clamp(t_opt_best.+rand(Uniform(minimum(ξ_0),maximum(ξ_0))),minimum(ξ_0),maximum(ξ_0)));
                result = Newton(combined_fun, combined_gradient, vars0; maxEval=max_iter, εAccuracy=ε)
                x_opt, t_opt = unpack_vars(result.xMin)
            end

            # Find the optimal q (Risk function) for the new optimal x only if the simplex and stochastic dominance constraints are satisfied.
            if convergence_simplex(x_opt) + convergence_SD(x_opt, t_opt, ξ, ξ_0, p, p_ξ, p_ξ_0) < ε
            q_opt = optimize_min_risk(vec(x_opt'*ξ),p_ξ;β=β,r=r,ε=1e-5,max_iter=50,n_particles=1000,max_ipot=50)
            end
            
            # Stopping condition:-
                # First: Check for at least one feasible point; if none is found, dominance is not possible  
                # Second: If the first condition is satisfied and a better objective is found, update the optimal values
                # Third: Try up to max_ipot times to look for further improvement
            if  ConvergenceRate_one_feasible_point > ε
                println("No single set of portfolio allocations exists that satisfies Stochastic Dominance of order $SDorder")
                    if verbose 
                        println("Consider increasing max_iter and n_particles, or reducing ε (current residuals: $ConvergenceRate_one_feasible_point), or trying for $SDorder+1.")
                        println("Simplex Constraints residuals: ", convergence_simplex(x_opt_best))
                        println("Stochastic Dominance Constraints residuals: ", convergence_SD(x_opt_best, t_opt_best, ξ, ξ_0, p, p_ξ, p_ξ_0))    
                    end 
                break
            elseif riskfunction(q_opt, vec(x_opt'*ξ), r, p_ξ, β) < riskfunction(q_opt_best, vec(x_opt_best'*ξ), r, p_ξ, β) && convergence_simplex(x_opt) + convergence_SD(x_opt, t_opt, ξ, ξ_0, p, p_ξ, p_ξ_0) < ε
    
                x_opt_best, q_opt_best, t_opt_best = x_opt, q_opt, t_opt
                failure_count = 0  # Reset failure count when the condition is met
    
            else
                failure_count += 1  # Increment failure count when the condition is not met
                if failure_count >= max_ipot
                    if verbose 
                        println("Simplex Constraints residuals: ", convergence_simplex(x_opt_best))
                        println("Stochastic Dominance Constraints residuals: ", convergence_SD(x_opt_best, t_opt_best, ξ, ξ_0, p, p_ξ, p_ξ_0))    
                    end 
                    if plot
                        B_q = optimize_min_risk(ξ_0,p_ξ_0;β=β,r=r,ε=1e-5,max_iter=50,n_particles=1000,max_ipot=50)
                        plot_min_return(x_opt_best,q_opt_best,B_q,ξ, p_ξ, ξ_0, p_ξ_0,p,r,β)
                    end
                    break  # Break the outer while loop after max_ipot consecutive failures                
                end
            end
    end
   return x_opt_best, q_opt_best, t_opt_best     
end

"""
    optimize_min_risk(ξ_0::AbstractVector{T}, p_ξ_0::AbstractVector{T}; 
        β::T=0.5, r::T=2, ε::Float64=1e-7, max_iter::Integer=50, 
        n_particles::Integer=100, max_ipot::Integer=50) where {T<:Number}

Performs **risk minimization** using Particle Swarm Optimization (PSO) to find the optimal risk threshold parameter `q`.

# Arguments
- `ξ_0::AbstractVector{<:Number}`: Benchmark return scenarios.
- `p_ξ_0::AbstractVector{<:Number}`: Probability distribution over benchmark scenarios.

# Keyword Arguments
- `β::Number = 0.5`: Confidence level for risk measure.
- `r::Number = 2.0`: Norm of the risk measure.
- `ε::Float64 = 1e-7`: Convergence tolerance.
- `max_iter::Integer = 50`: Maximum number of iterations for optimization.
- `n_particles::Integer = 100`: Number of particles for PSO.
- `max_ipot::Integer = 50`: Maximum number of attempts for interior point optimization.

# Returns
- `q_opt_best::Float64`: The optimized risk threshold parameter.

# Example
```julia
# Define benchmark return scenarios
ξ_0 = [0.01, 0.04, 0.00]  # Benchmark returns in different scenarios

# Define probability distribution for benchmark scenarios
p_ξ_0 = [0.4, 0.4, 0.2]   # Probabilities for benchmark scenarios

# Set risk parameters
β = 0.5
r = 2.0

# Solve the risk minimization problem
q_opt_best = optimize_min_risk(ξ_0, p_ξ_0; β=β, r=r, ε=1e-6, max_iter=100)
```
"""
function optimize_min_risk(ξ_0::AbstractVector{T}, 
    p_ξ_0::AbstractVector{T}; 
    β::T=0.5, 
    r::T=2, 
    ε::Float64=1e-7, 
    max_iter::Integer=50, 
    n_particles::Integer=100, 
    max_ipot::Integer=50) where {T<:Number}

    
    # PSO Initialization
    q=rand()

    max_ipot = max_ipot  # Maximum times IP Optimization tries
    
    # Objective function (risk measure) for optimization 
    function combined_fun(vars)
      q = unpack_vars(vars)
      result= riskfunction(q, ξ_0, r, p_ξ_0, β)
        return result
    end
    # Function to unpack variables from optimization vector
    function unpack_vars(vars)
        q = vars[1]              # Next element is q
        return q
    end
  
    # Derivative of the objective function (risk measure) stopping condition 
    function grad_lagrangian_q_risk(q, ξ_0, r, p_ξ_0, β)
        term1 =(safe_exponent.(max.(-ξ_0 .- q, 0), r)'* p_ξ_0).^((1/r) - 1)
        if term1 == Inf 
            term1 = 0
        end
        mainterm = 1 .- (1 / (1 - β)) * term1 *
                    (safe_exponent.(max.(-ξ_0 .- q, 0), r-1)' * p_ξ_0)
        return norm(mainterm)
    end
    # Convert the derivative of the objective function (risk measure) into a single-input function (q_opt) for optimization
    function combined_Lagfun(vars)
        q = unpack_vars(vars)
        result= grad_lagrangian_q_risk(q, ξ_0, r, p_ξ_0, β)
        return result
      end
    
    # Initialize optimization variables
    vars0 = vcat(q)
    # Optimize using PSO (initial) 
    result= PSO(combined_fun,vars0;n_particles = n_particles,maxEval = max_iter, stoppingcondition=true,stoppingconditionfun=combined_Lagfun)
    # Result after first try 
    q_opt_best = unpack_vars(result.xMin)

    failure_count = 0   # Counter for tracking failures
    # Try up to max_ipot times or until the derivative of the objective function converges to zero using PSO to find the best optimal point q
    while true
        vars0 = vcat(q_opt_best)
        result= PSO(combined_fun,vars0;n_particles = n_particles,maxEval = max_iter,Initialization=true,stoppingcondition=true,stoppingconditionfun=combined_Lagfun)
        q_opt= unpack_vars(result.xMin)
        if grad_lagrangian_q_risk(q_opt, ξ_0, r, p_ξ_0, β) < grad_lagrangian_q_risk(q_opt_best, ξ_0, r, p_ξ_0, β)
            q_opt_best= q_opt
        else 
            failure_count += 1 
            if grad_lagrangian_q_risk(q_opt, ξ_0, r, p_ξ_0, β) < ε 
                break 
            elseif failure_count >= max_ipot
                break 
            end            
        end
    end    
    return  q_opt_best
end



