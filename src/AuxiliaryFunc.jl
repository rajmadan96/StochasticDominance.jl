module AuxiliaryFunc

using LinearAlgebra, Plots 

export safe_exponent, simplex1, simplex2, g_p_minus_1, g_p, g_bar, boundary_moments, CheckConvergenceSimplex, CheckConvergenceSD, MeanReturn, plotOptimalAssetAllocationMeanReturn, RiskFunction, BenchmarkRiskFunction, plotOptimalAssetAllocationRiskFunction

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
- `x::AbstractVector{<:Number}`: Portfolio weights (simplex).
- `ξ::AbstractVector{<:Number}`: Portfolio returns.
- `p_ξ::AbstractVector{<:Number}`: Probability of portfolio `ξ`.
- `p_ξ_0::AbstractVector{<:Number}`: Probability of benchmark `ξ_0`.
- `ξ_0::AbstractVector{<:Number}`: Benchmark returns.
- `p::Integer`: The highest moment order.

# Returns
- A vector of length `p` containing computed moment constraints:
  - For `k = 1`, compares the mean of the portfolio and benchmark.
  - For `k ≥ 2`, computes and compares the central moments.
  - If the portfolio moment does not exceed the benchmark moment, the value is set to zero.
  - Otherwise, the result stores the norm of the difference.

# Examples
```julia
boundary_moments([0.7, 0.3], [2  3; 5 6], [0.5,0.5], [0.2,0.3,0.5], [2,3,4], 3)  
```
"""
function boundary_moments(x, ξ, p_ξ,p_ξ_0,ξ_0, p)
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
            centered_vec = vec_xξ .- mean_xξ
            lhs = dot(centered_vec .^ k, p_ξ)
            rhs = p_ξ_0'*(max.((maximum(ξ_0) .- ξ_0) .^ k,0))
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

"""
    CheckConvergenceSimplex(x)

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
CheckConvergenceSimplex([0.4, 0.6])    # Expected output: 0.0 (valid simplex)
CheckConvergenceSimplex([0.5, 0.7])    # Nonzero output (sum > 1)
CheckConvergenceSimplex([0.5, -0.2])   # Nonzero output (contains negative values)
```
"""
function CheckConvergenceSimplex(x)
    simplex_constraints = vcat(
        simplex1(x),
        simplex2(x)
    )
    return norm(simplex_constraints)
end

"""
    CheckConvergenceSD(x, t, ξ, ξ_0, p, p_ξ, p_ξ_0)

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
CheckConvergenceSD([0.7, 0.3], 5, [2  3; 5 6], [2,3,4], 3, [0.5,0.5], [0.2,0.3,0.5])  
# Computes the norm of stochastic dominance constraint violations.
```
"""
function CheckConvergenceSD(x, t, ξ, ξ_0, p, p_ξ, p_ξ_0;)
    stochastic_dominance_constraints = vcat(
        g_p(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0),
        g_p_minus_1(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0),
        g_bar(x, ξ, ξ_0, p, p_ξ, p_ξ_0),
        boundary_moments(x, ξ, p_ξ,p_ξ_0,ξ_0, p)
    )
    return norm(stochastic_dominance_constraints) 
end

"""
    MeanReturn(x, ξ, p_ξ)

Compute the expected return of a portfolio given asset returns and probabilities.

# Arguments
- `x::AbstractVector{<:Number}`: Portfolio weights (simplex).
- `ξ::AbstractMatrix{<:Number}`: Matrix of asset returns (each column represents a scenario).
- `p_ξ::AbstractVector{<:Number}`: Probability distribution over the scenarios.

# Returns
- A numeric value representing the expected portfolio return 

# Examples
```julia
MeanReturn([0.7, 0.3], [2  3; 5 6], [0.5,0.5])  
# Computes the expected portfolio return.
```
"""
function MeanReturn(x,ξ,p_ξ)
   return  dot(vec((p_ξ)'*ξ'),x)
end

"""
    plotOptimalAssetAllocationMeanReturn(x, ξ, p_ξ, ξ_0, p_ξ_0)

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
plotOptimalAssetAllocationMeanReturn([0.4, 0.6], [2  3; 5 6], [0.5,0.5], [2,3], [0.2,0.8])  
# Generates a pie chart of asset allocation and displays mean return vs. benchmark.
```
"""
function plotOptimalAssetAllocationMeanReturn(x, ξ, p_ξ, ξ_0, p_ξ_0)
    meanreturn = round(MeanReturn(x, ξ, p_ξ), digits=3)
    benchmark = round(p_ξ_0' * ξ_0, digits=3)
    x_percentage= round.(abs.(x)./sum(abs.(x))*100,digits=1)
    labels = ["Asset $i: $(x_percentage[i])%" for i in 1:length(x)]
    AssetPlot=pie(labels, x, title="Optimal Asset Allocation (SD Order $(p+1))", legend=legend=(0.01, 0.85))
    annotate!(0, -1.1, text("Mean Return: $meanreturn; Benchmark: $benchmark", :center, 10))
    
    display(plot(AssetPlot))
end

"""
    RiskFunction(x, q, ξ, p, p_ξ, β)

Compute the risk function based on a distortion risk measure, incorporating 
higher-order moments and tail risk.

# Arguments
- `x::AbstractVector{<:Number}`: Portfolio weights.
- `q::Number`: Risk threshold parameter.
- `ξ::AbstractMatrix{<:Number}`: Matrix of asset returns (each column represents a scenario).
- `p::Number`: Exponent parameter controlling tail risk sensitivity.
- `p_ξ::AbstractVector{<:Number}`: Probability distribution over the scenarios.
- `β::Number`: Confidence level (closer to 1 indicates higher aversion).

# Returns
- A numeric value representing the risk function  
# Examples
```julia
RiskFunction([0.4, 0.6], 1.5, [2  3; 5 6], 3, [0.5,0.5], 0.95)  
# Computes the risk function for the given parameters.
```
"""
function RiskFunction(x, q, ξ, p, p_ξ, β)
    return q + (1 / (1 - β)) * (sum(p_ξ[i] * safe_exponent(max(-dot(x, ξ[:, i]) - q, 0), p) for i in 1:size(ξ, 2)))^(1/p)
end

"""
    BenchmarkRiskFunction(q, ξ_0, p, p_ξ_0, β)

Compute the risk function for the benchmark portfolio using a distortion risk measure, 
capturing higher-order moments and tail risk.

# Arguments
- `q::Number`: Risk threshold parameter.
- `ξ_0::AbstractVector{<:Number}`: Benchmark returns.
- `p::Number`: Exponent parameter controlling tail risk sensitivity.
- `p_ξ_0::AbstractVector{<:Number}`: Probability distribution over the benchmark scenarios.
- `β::Number`: Confidence level (closer to 1 indicates higher aversion).

# Returns
- A numeric value representing the benchmark risk function
# Examples
```julia
BenchmarkRiskFunction(1.5, [2, 3, 4], 3, [0.2, 0.3, 0.5], 0.95)  
# Computes the benchmark risk function for the given parameters.
```
"""
function BenchmarkRiskFunction(q, ξ_0, p, p_ξ_0, β)
    return q + (1 / (1 - β)) * (sum(p_ξ_0[i] * safe_exponent(max(-ξ_0[i] - q, 0), p) for i in 1:length(ξ_0)))^(1/p)
end

"""
    plotOptimalAssetAllocationRiskFunction(x, q, ξ, p_ξ, ξ_0, p_ξ_0, p)

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
function plotOptimalAssetAllocationRiskFunction(x,q,ξ, p_ξ, ξ_0, p_ξ_0,p)
    OptRisk  = round( RiskFunction(x, q, ξ, p, p_ξ, β), digits=3)
    BenchmarkRisk = round(BenchmarkRiskFunction(q, ξ_0, p, p_ξ_0, β), digits=3)
    x_percentage= round.(abs.(x)./sum(abs.(x))*100,digits=1)
    labels = ["Asset $i: $(x_percentage[i])%" for i in 1:length(x)]
    
    AssetPlot=pie(labels, x, title="Optimal Asset Allocation (SD Order $(p+1))", legend=legend=(0.01, 0.85))
    annotate!(0, -1.1, text("Optimal Risk: $OptRisk; Benchmark Risk: $BenchmarkRisk", :center, 10))
    display(plot(AssetPlot))
end

end # end of the module