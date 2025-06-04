using LinearAlgebra


"""
    verify_dominance(Y, X,SDorder; p_X, p_Y,ε=1e-5,verbose=false)

Verify if `Y` stochastically dominates `X` at a given stochastic oder 

# Arguments
- `Y::AbstractVector{<:Number}`: Portfolio outcomes (returns).
- `X::AbstractVector{<:Number}`: Benchmark outcomes.
- `p_X::AbstractVector{<:Number}`: Probabilities associated with `X`.
- `p_Y::AbstractVector{<:Number}`: Probabilities associated with `Y`.
- `SDorder::Integer`: Order of stochastic dominance (SD).

# Keyword Arguments
- `ε::Float64=1e-5`: Convergence tolerance for the dominance check.
- `verbose::Bool=false`: If `true`, print detailed accuracy or residuals.

# Returns
- A message indicating whether `Y` stochastically dominates `X` at order `SDorder`.
- `true` if `Y` stochastically dominates `X` at the given order.
- `false` otherwise.

# Examples
```julia
Y = [3, 5, 7]
X = [2, 4, 6]
p_Y = [0.3, 0.4, 0.3]
p_X = [0.2, 0.5, 0.3]
SDorder = 2

verify_dominance(Y, X, SDorder; p_X, p_Y, ε=1e-8,verbose=true)
```
"""
function verify_dominance(Y::AbstractVector{T}, 
    X::AbstractVector{T}, 
    SDorder::S; 
    p_Y=nothing, 
    p_X=nothing, 
    ε::Float64=1e-5, 
    verbose::Bool=false) where {T<:Number, S<:Number}
    # Validate inputs
    if SDorder == 1
        error("First order is not available. Try order: $(SDorder+1)")
    elseif length(Y) ≠ length(X)
        error("The length of portfolio and benchmark scenarios are not equal. They must be equal to make the stochastic dominance comparison.")
    end

    # Handling probabilities of X and Y
    if p_X !== nothing
        if length(p_X) == length(X)
            X_rand = randomVariable(X, p_X)
        else
            error("The length of benchmark scenarios and their respective probabilities are not equal. They must be equal.")
        end
    else
        X_rand = randomVariable(X)
    end
    p_X= X_rand.probabilities
    if p_Y !== nothing
        if length(p_Y) == length(Y)
            Y_rand = randomVariable(Y, p_Y)
        else
            error("The length of portfolio scenarios and their respective probabilities are not equal. They must be equal.")
        end
    else
        Y_rand = randomVariable(Y)
    end
    p_Y= Y_rand.probabilities

    # Define variables
    ξ = I(length(Y))  # Identity matrix of appropriate size
    x=Y
    ξ_0 = X
    ξ = I(length(Y))  # Identity matrix of appropriate size
    p_ξ = p_Y
    p_ξ_0 = p_X
    p=SDorder-1.0
    
    # Compute the conditions to check whether Y dominates X or not
    stochastic_dominance_constraints = vcat(
            g_bar(x, ξ, ξ_0, p, p_ξ, p_ξ_0),
        boundary_moments(x, ξ, p_ξ,p_ξ_0,ξ_0, p))

    # Output    
    if norm(stochastic_dominance_constraints)  ≤ ε
        if verbose
            accuracy = (1-norm(stochastic_dominance_constraints))*100
            println("Y dominates X in stochastic order $SDorder with $accuracy % accuracy")
        end
        println("Y dominates X in stochastic order $SDorder")
        return true 
    else
        if verbose
                residual = norm(stochastic_dominance_constraints)
                println("Y doesn't dominates X in stochastic order $SDorder with residual $residual") 
        end
        println("Y doesn't dominates X in stochastic order $SDorder")
        return false
    end
end






