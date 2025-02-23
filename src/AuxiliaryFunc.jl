module AuxiliaryFunc

using LinearAlgebra, Plots 

export safe_exponent, simplex1, simplex2, g_p_minus_1, g_p, g_bar, boundary_moments, 
       convergence_simplex, convergence_SD, expected_portfolio_return, 
       plot_max_return, riskfunction_asset_allocation, riskfunction, 
       plot_min_return, randomVariable

"""
    safe_exponent(base, exponent)

Compute `base^exponent` while handling the indeterminate case `0^0`.

# Arguments
- `base::Number`: The base value.
- `exponent::Number`: The exponent value.

# Returns
- `base^exponent` if `(base, exponent) ≠ (0, 0)`.
- `0.0` if `base == 0` and `exponent == 0`.

# Examples
```julia
safe_exponent(2, 3)   # 8
safe_exponent(5, 0)   # 1
safe_exponent(0, 5)   # 0
safe_exponent(0, 0)   # 0.0
safe_exponent(2, -2)  # 0.25
```
"""
function safe_exponent(base, exponent)
    if base == 0 && exponent == 0
        return 0.0
    else
        return base^exponent
    end
end

"""
    simplex1(x)

Compute `1 - sum(x)`, which represents the residual of the sum of elements in `x` from 1.

# Arguments
- `x::AbstractArray{<:Number}`: An array of numerical values.

# Returns
- A numeric value representing `1 - sum(x)`.

# Examples
```julia
simplex1([0.2, 0.3, 0.1])  # 0.4
simplex1([0.5, 0.5])       # 0.0
simplex1([1.2, -0.3])      # 0.1
```
"""
function simplex1(x)
    result = 1 - sum(x)
    return result
end

"""
    simplex2(x)

Compute the element-wise product of `x` and a boolean mask indicating negative values in `x`. 
This function extracts the negative components of `x`.

# Arguments
- `x::AbstractArray{<:Number}`: An array of numerical values.

# Returns
- An array where negative elements of `x` are preserved, and non-negative elements are replaced with `0`.

# Examples
```julia
simplex2([1, -2, 3, -4])   # [0, -2, 0, -4]
simplex2([0.5, -0.1, 0])   # [0.0, -0.1, 0.0]
simplex2([-3, -2, -1])     # [-3, -2, -1]
```
"""
function simplex2(x)
    result = x .* (x .< 0)
    return x .* (x .< 0)
end

"""
    g_p_minus_1(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0)

Compute the derivative of Higher order stochastic domiance constraint for single t, where t ∈ R. 

# Arguments
- `t::Number`: A scalar threshold value.
- `x::AbstractVector{<:Number}`: portfolio weights (simplex).
- `ξ::AbstractVector{<:Number}`: portfolio.
- `ξ_0::Number`: Benchmark returns.
- `p::Number`: The exponent parameter.
- `p_ξ::Number`: Probability of portfolio `ξ`.
- `p_ξ_0::Number`: Probability of benchmark  `ξ_0`.

# Returns
- A numeric value computed as `term1[1] - term2[1]`, where:
  - `term1` involves `safe_exponent` applied to `max(t - x' * ξ, 0)`, scaled by `p_ξ`.
  - `term2` involves `safe_exponent` applied to `max(t - ξ_0, 0)`, scaled by `p_ξ_0`.

# Examples
```julia
g_p_minus_1(5, [0.7, 0.3], [2  3; 5 6], [2,3], 3, [0.5,0.5], [0.2,0.8])  
```
"""
function g_p_minus_1(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0)
    term1 = safe_exponent.(max.(t .- x' * ξ, 0), p-1) * p_ξ
    term2 = safe_exponent.(max.(t .- ξ_0, 0), p-1)' * p_ξ_0
    return term1[1] - term2[1]
end


"""
    g_p(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0)

Compute the higher-order stochastic dominance constraint for a single `t`, where `t ∈ R`.

# Arguments
- `t::Number`: A scalar threshold value.
- `x::AbstractVector{<:Number}`: Portfolio weights (simplex).
- `ξ::AbstractVector{<:Number}`: Portfolio returns.
- `ξ_0::Number`: Benchmark returns.
- `p::Number`: The exponent parameter.
- `p_ξ::Number`: Probability of portfolio `ξ`.
- `p_ξ_0::Number`: Probability of benchmark `ξ_0`.

# Returns
- A numeric value computed as `term1[1] - term2[1]`, where:
  - `term1` involves `safe_exponent` applied to `max(t - x' * ξ, 0)`, scaled by `p_ξ`.
  - `term2` involves `safe_exponent` applied to `max(t - ξ_0, 0)`, scaled by `p_ξ_0`.

# Examples
```julia
g_p(5, [0.7, 0.3], [2  3; 5 6], [2,3], 3, [0.5,0.5], [0.2,0.8])  
```
"""
function g_p(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0)
    term1 = (safe_exponent.(max.(t .- x' * ξ, 0), p) * p_ξ)
    term2 = (safe_exponent.(max.(t .- ξ_0, 0), p)' * p_ξ_0)
    return term1[1] - term2[1]
end

"""
    g_bar(x, ξ, ξ_0, p, p_ξ, p_ξ_0)

Compute the verification of the higher-order stochastic dominance constraint by evaluating `g_p` over sorted benchmark values.

# Arguments
- `x::AbstractVector{<:Number}`: Portfolio weights (simplex).
- `ξ::AbstractVector{<:Number}`: Portfolio returns.
- `ξ_0::AbstractVector{<:Number}`: Benchmark returns.
- `p::Number`: The exponent parameter.
- `p_ξ::AbstractVector{<:Number}`: Probability of portfolio `ξ`.
- `p_ξ_0::AbstractVector{<:Number}`: Probability of benchmark `ξ_0`.

# Returns
- A numeric value computed as `max(max_value, 0.0)`, where:
  - `max_value` is the maximum of `g_p(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0)` over sorted `ξ_0`.
  - If `max_value` is negative, the function returns `0.0`.
  - If `max_value` is positive, the selected portfolio or the portfolio with assigned weights does not satisfy the dominance condition."

# Examples
```julia
g_bar([0.7, 0.3], [2  3; 5 6], [2,3,4], 3, [0.5,0.5], [0.2,0.3,0.5])  
```
"""
function g_bar(x, ξ, ξ_0, p, p_ξ, p_ξ_0)
    # Sort ξ_0
    sorted_ξ_0 = sort(ξ_0)

    # Collect results into an array directly
    results = [
        g_p(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0) for t in sorted_ξ_0
    ]
    
    # Calculate the maximum value
    max_value = maximum(results) #maximum(results)
    
    # Return the exact maximum if non-negative, otherwise return zero
    return max(max_value, 0.0)
end

"""
boundary_moments(x, ξ, p_ξ, p_ξ_0, ξ_0, p)

Compute the boundary moments for higher-order stochastic dominance constraints by comparing 
the moments of the portfolio and benchmark distributions.

# Arguments
- `x::AbstractVector{<:Number}`: Portfolio weights (simplex constraint).
- `ξ::AbstractMatrix{<:Number}`: Portfolio returns.
- `p_ξ::AbstractVector{<:Number}`: Probability distribution of the portfolio returns `ξ`.
- `p_ξ_0::AbstractVector{<:Number}`: Probability distribution of the benchmark returns `ξ_0`.
- `ξ_0::AbstractVector{<:Number}`: Benchmark returns.
- `p::Real`: The highest moment order. If `p` is an integer, `boundary_moments_integer` is used. If `p` is non-integer, `boundary_moments_noninteger` is used.

# Returns
- A vector containing computed moment constraints:
- For `k = 1`: Compares the mean of the portfolio and benchmark.
- For `k ≥ 2`: Computes and compares the central moments.
- If the portfolio moment does not exceed the benchmark moment, the value is set to zero.
- Otherwise, the result stores the norm of the difference.

# Methodology
- Integer Case (`p` is an integer):
- Computes moments from order 1 order `p`.
- Uses the function `boundary_moments_integer`.
- Non-Integer Case (`p` is a non-integer):
- Computes moments from order (0,1) order `p`
- Uses the function `boundary_moments_noninteger`.

# Examples
```julia
# Example 1: Integer order moments
boundary_moments([0.7, 0.3], [2 3; 5 6], [0.5, 0.5], [0.2, 0.3, 0.5], [2, 3, 4], 3)

# Example 2: Non-integer order moments
boundary_moments([0.7, 0.3], [2 3; 5 6], [0.5, 0.5], [0.2, 0.3, 0.5], [2, 3, 4], 2.7)
```
"""
function boundary_moments(x, ξ, p_ξ, p_ξ_0, ξ_0, p)
    if isinteger(p)
        return boundary_moments_integer(x, ξ, p_ξ, p_ξ_0, ξ_0, p)
    else
        return boundary_moments_noninteger(x, ξ, p_ξ, p_ξ_0, ξ_0, p)
    end
end

function boundary_moments_integer(x, ξ, p_ξ,p_ξ_0,ξ_0, p)
    # Convert p to an integer, pre-allocate result vector
    p = Int(p)
    
    # Determine the type for initialization
    result_type = promote_type(eltype(x' * ξ), eltype(p_ξ), eltype(ξ_0))
    result = zeros(result_type, p)  # Pre-allocate the result vector

    # Compute the mean of the random variable vec(x' * ξ)
    vec_xξ = vec(x' * ξ)
    mean_xξ = p_ξ'*vec_xξ  # Mean of vec(x' * ξ)

    # Compute the mean-centered random variable vec(x' * ξ) 
    centered_vec = max.(maximum(ξ_0) .- vec_xξ , 0)

    # Compute results for each k = 1, 2, ..., p
    for k in 1:p
        if k == 1
            # For k = 1, compare the means
            lhs = -mean_xξ
            rhs = -p_ξ_0'*ξ_0
        else
            # For k >= 2, compute central moments
            lhs = dot(centered_vec .^ k, p_ξ)
            rhs = p_ξ_0'*((maximum(ξ_0) .- ξ_0) .^ k)
        end

        # Compare lhs and rhs
        if lhs ≤ rhs
            result[k] = zero(result_type)
        else
            result[k] = norm(lhs - rhs)
        end
    end
    return result
end

function boundary_moments_noninteger(x, ξ, p_ξ,p_ξ_0,ξ_0, p)
    # Convert p to an integer, pre-allocate result vector
    
    temp_p = Int(floor(p))
    
    # Determine the type for initialization
    result_type = promote_type(eltype(x' * ξ), eltype(p_ξ), eltype(ξ_0))
    result = zeros(result_type, temp_p+1)  # Pre-allocate the result vector

    # Compute the mean of the random variable vec(x' * ξ)
    vec_xξ = vec(x' * ξ)
    mean_xξ = p_ξ'*vec_xξ  # Mean of vec(x' * ξ)

    # Compute the mean-centered random variable vec(x' * ξ) 
    centered_vec = max.(maximum(ξ_0) .- vec_xξ , 0)

    # Compute results for each k = 1, 2, ..., p
    for k in 0:temp_p
            lhs = dot(centered_vec .^(p-k), p_ξ)
            rhs = p_ξ_0'*((maximum(ξ_0) .- ξ_0) .^(p-k))
                # Compare lhs and rhs
        if lhs ≤ rhs
            result[k+1] = zero(result_type)
        else
            result[k+1] = norm(lhs - rhs)
        end
    end
    return result
end

"""
    convergence_simplex(x)

Evaluate the convergence of a given vector `x` to a simplex by computing the norm 
of its deviation from simplex constraints.

# Arguments
- `x::AbstractVector{<:Number}`: A vector representing portfolio weights or probabilities.

# Returns
- A numeric value representing the norm of simplex constraint violations:
  - Uses `simplex1(x)` to check if the sum of elements deviates from 1.
  - Uses `simplex2(x)` to check for negative elements.
  - Returns the norm of the concatenated constraint violations.

# Examples
```julia
convergence_simplex([0.4, 0.6])    # Expected output: 0.0 (valid simplex)
convergence_simplex([0.5, 0.7])    # Nonzero output (sum > 1)
convergence_simplex([0.5, -0.2])   # Nonzero output (contains negative values)
```
"""
function convergence_simplex(x)
    simplex_constraints = vcat(
        simplex1(x),
        simplex2(x)
    )
    return norm(simplex_constraints)
end

"""
    convergence_SD(x, t, ξ, ξ_0, p, p_ξ, p_ξ_0)

Evaluate the convergence of a given portfolio `x` to satisfy stochastic dominance constraints 
by computing the norm of constraint violations.

# Arguments
- `x::AbstractVector{<:Number}`: Portfolio weights (simplex).
- `t::Number`: A scalar threshold value for dominance constraints.
- `ξ::AbstractVector{<:Number}`: Portfolio returns.
- `ξ_0::AbstractVector{<:Number}`: Benchmark returns.
- `p::Integer`: The order of stochastic dominance.
- `p_ξ::AbstractVector{<:Number}`: Probability of portfolio `ξ`.
- `p_ξ_0::AbstractVector{<:Number}`: Probability of benchmark `ξ_0`.

# Returns
- A numeric value representing the norm of stochastic dominance constraint violations:
  - Uses `g_p(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0)` to evaluate the primary dominance constraint.
  - Uses `g_p_minus_1(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0)` to check lower-order dominance.
  - Uses `g_bar(x, ξ, ξ_0, p, p_ξ, p_ξ_0)` to compute the upper bound constraint.
  - Uses `boundary_moments(x, ξ, p_ξ, p_ξ_0, ξ_0, p)` to ensure moment-based dominance.
  - Returns the norm of these combined constraints.

# Examples
```julia
convergence_SD([0.7, 0.3], 5, [2  3; 5 6], [2,3,4], 3, [0.5,0.5], [0.2,0.3,0.5])  
# Computes the norm of stochastic dominance constraint violations.
```
"""
function convergence_SD(x, t, ξ, ξ_0, p, p_ξ, p_ξ_0;)
    stochastic_dominance_constraints = vcat(
        g_p(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0),
        g_p_minus_1(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0),
        g_bar(x, ξ, ξ_0, p, p_ξ, p_ξ_0),
        boundary_moments(x, ξ, p_ξ,p_ξ_0,ξ_0, p)
    )
    return norm(stochastic_dominance_constraints) 
end

"""
    expected_portfolio_return(x, ξ, p_ξ)

Compute the expected return of a portfolio given asset returns and probabilities.

# Arguments
- `x::AbstractVector{<:Number}`: Portfolio weights (simplex).
- `ξ::AbstractMatrix{<:Number}`: Matrix of asset returns (each column represents a scenario).
- `p_ξ::AbstractVector{<:Number}`: Probability distribution over the scenarios.

# Returns
- A numeric value representing the expected portfolio return 

# Examples
```julia
expected_portfolio_return([0.7, 0.3], [2  3; 5 6], [0.5,0.5])  
# Computes the expected portfolio return.
```
"""
function expected_portfolio_return(x,ξ,p_ξ)
   return  dot(vec((p_ξ)'*ξ'),x)
end

"""
    plot_max_return(x, ξ, p_ξ, ξ_0, p_ξ_0)

Plot the optimal asset allocation as a pie chart, displaying portfolio weights 
and comparing mean return to a benchmark.

# Arguments
- `x::AbstractVector{<:Number}`: Portfolio weights (simplex).
- `ξ::AbstractMatrix{<:Number}`: Matrix of asset returns (each column represents a scenario).
- `p_ξ::AbstractVector{<:Number}`: Probability distribution over the portfolio scenarios.
- `ξ_0::AbstractVector{<:Number}`: Benchmark returns.
- `p_ξ_0::AbstractVector{<:Number}`: Probability distribution over the benchmark scenarios.

# Returns
- A pie chart illustrating the portfolio weights.
- Annotates the mean return and benchmark return.

# Examples
```julia
plot_max_return([0.4, 0.6], [2  3; 5 6], [0.5,0.5], [2,3], [0.2,0.8])  
# Generates a pie chart of asset allocation and displays mean return vs. benchmark.
```
"""
function plot_max_return(x, ξ, p_ξ, ξ_0, p_ξ_0)
    portfolio_return = round(expected_portfolio_return(x, ξ, p_ξ), digits=3)
    benchmark = round(p_ξ_0' * ξ_0, digits=3)
    x_percentage= round.(abs.(x)./sum(abs.(x))*100,digits=1)
    labels = ["asset $i: $(x_percentage[i])%" for i in 1:length(x)]
    AssetPlot=pie(labels, x, title="Optimal asset allocation (SD order $(p+1))", legend=legend=(0.01, 0.85))
    annotate!(0, -1.1, text("maximized expected return: $portfolio_return%\n return of benchmark: $benchmark%", :center, 10, :black))
    display(plot(AssetPlot))
end


"""
    riskfunction_asset_allocation(x, q, ξ, r, p_ξ, β)

Calculates the risk-adjusted portfolio return based on quantile-based risk measures.

# Arguments
- `x`: Portfolio weights (decision variables to optimize).
- `q`: Quantile threshold for the risk measure.
- `ξ`: Portfolio scenarios (matrix of asset returns).
- `r`: Norm parameter.
- `p_ξ`: Probability weights for the scenarios in `ξ`.
- `β`: Risk aversion parameter (closer to 1 indicates higher aversion).

# Returns
- The risk-adjusted return, which accounts for portfolio performance and risk under extreme loss scenarios.
"""
function riskfunction_asset_allocation(x, q, ξ, r, p_ξ, β)
    return q + (1 / (1 - β)) * (sum(p_ξ[i] * safe_exponent(max(-dot(x, ξ[:, i]) - q, 0),r) for i in 1:size(ξ, 2)))^(1/r)
end

"""
    riskfunction(q, ξ_0, r, p_ξ_0, β)

Calculates the risk-adjusted return based on quantile-based risk measures.

# Arguments
- `q`: Quantile threshold for the risk measure.
- `ξ_0`: Random variable 
- `r`: Norm parameter.
- `p_ξ_0`: Probability weights for the scenarios in `ξ`.
- `β`: Risk aversion parameter (closer to 1 indicates higher aversion).

# Returns
- The risk-adjusted return, risk under extreme loss scenarios.
"""
function riskfunction(q, ξ_0, r, p_ξ_0, β)
    return q + (1 / (1 - β)) * ((sum(p_ξ_0[i] * safe_exponent(max.(-ξ_0[i] .- q, 0), r) for i in 1:length(ξ_0)))[1])^(1/r)
end


"""
    plot_min_return(x, q, ξ, p_ξ, ξ_0, p_ξ_0, p)

Plot the optimal asset allocation as a pie chart, displaying portfolio weights and comparing 
the portfolio risk function to the benchmark risk function.

# Arguments
- `x::AbstractVector{<:Number}`: Portfolio weights (simplex).
- `q::Number`: Risk threshold parameter.
- `ξ::AbstractMatrix{<:Number}`: Matrix of asset returns (each column represents a scenario).
- `p_ξ::AbstractVector{<:Number}`: Probability distribution over the portfolio scenarios.
- `ξ_0::AbstractVector{<:Number}`: Benchmark returns.
- `p_ξ_0::AbstractVector{<:Number}`: Probability distribution over the benchmark scenarios.
- `p::Number`: Exponent parameter controlling tail risk sensitivity.

# Returns
- A pie chart illustrating the portfolio weights.
- Annotates the optimal portfolio risk and benchmark risk.

# Examples
```julia
plotOptimalAssetAllocationRiskFunction([0.4, 0.6], 1.5, [2  3; 5 6], [0.5,0.5], [2,3], [0.2,0.8], 3)  
```
"""
function plot_min_return(x,q,B_q,ξ, p_ξ, ξ_0, p_ξ_0,p,r)
    OptRisk  = round(riskfunction_asset_allocation(x, q, ξ, r, p_ξ, β), digits=3)
    BenchmarkRisk = round(riskfunction(B_q, ξ_0, r, p_ξ_0, β), digits=3)
    x_percentage= round.(abs.(x)./sum(abs.(x))*100,digits=1)
    labels = ["asset $i: $(x_percentage[i])%" for i in 1:length(x)]
    
    AssetPlot=pie(labels, x, title="Optimal asset allocation (SD order $(p+1))", legend=legend=(0.01, 0.85))
    annotate!(0, -1.1, text("optimal coherent risk-return →  portfolio: $OptRisk%, benchmark: $BenchmarkRisk%", :center, 10))
    #annotate!(0, -1.1, text("Optimal Risk: $OptRisk", :center, 10))
    display(plot(AssetPlot))
end

"""
    ontosimplex!(probabilities::Vector{Float64})

Projects the given probability vector onto the probability simplex, ensuring all elements are 
non-negative and sum to 1.

# Arguments
- `probabilities::Vector{Float64}`: A vector of probabilities that may not initially sum to 1 or contain negative values.

# Behavior
- If any probability is negative, it is set to 0.
- If the sum of probabilities is non-positive, a uniform distribution is assigned.
- If the probabilities do not sum to 1, they are normalized, and a random index is adjusted to ensure the sum equals 1.

# Returns
- Modifies `probabilities` in place to ensure it is a valid probability distribution.

# Examples
```julia
probabilities = [0.2, -0.1, 0.5, 0.3]
ontosimplex!(probabilities)
println(probabilities)  # Output: A valid probability vector summing to 1
"""

function ontosimplex!(probabilities::Vector{Float64})   # projects the probabilities onto the simplex
	for i= 1:length(probabilities)
		if probabilities[i] < 0.0 
			probabilities[i] = 0.0		# have them positive
		end
	end
	summ= sum(probabilities)
	if summ <= 0.0
		probabilities= fill!(probabilities, 1.0/ length(probabilities))
		summ= sum(probabilities)
		#@info "ontoSimplex: created probabilities." maxlog=1
	end
	if summ != 1
		probabilities./= summ
		tmpi= rand(1:length(probabilities))		# modify a random index
		summ= -sum(probabilities[1:tmpi-1])+ 1.0- sum(probabilities[tmpi+1:end])
		probabilities[tmpi]= max(0, summ)
		#@info "ontoSimplex: modified probabilities." maxlog=1
	end
	nothing		# modifies probabilities
end

struct randomVariable{T}
    states::Vector{T}              
    probabilities::Vector{Float64}
    length::Int                    

    function randomVariable(states::Vector{T}, probabilities::Vector{Float64} = fill(1.0 / length(states), length(states))) where T
        @assert length(states) == length(probabilities) "randomVariable: length of states and probabilities do not match."
        if any(probabilities .< 0.0)
            ontosimplex!(probabilities)
        elseif !(sum(probabilities) ≈ 1.0)
            ontosimplex!(probabilities)
        end
        new{T}(states, probabilities, length(states))
    end
end


end # end of the module
