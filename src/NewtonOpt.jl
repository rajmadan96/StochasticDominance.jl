module NewtonOpt

using LinearAlgebra, Random
Random.seed!(42)

export Newton

"""
    Newton(fun, funD, x0; maxEval=1000, εAccuracy=1e-7)

Perform Newton's method for root-finding or optimization.

# Arguments
- `fun::Function`: The function whose root or minimum is being sought.
- `funD::Function`: The derivative (Jacobian) of `fun`.
- `x0::Vector{Float64}`: Initial guess for the root or minimum.
- `maxEval::Int` (optional, default = `1000`): Maximum number of function evaluations.
- `εAccuracy::Float64` (optional, default = `1e-7`): Convergence tolerance.

# Returns
- A named tuple containing:
  - `xMin`: The estimated root or minimum point.
  - `fMin`: The function value at `xMin`.
  - `evalCount`: The number of function evaluations.

# Examples
```julia
f(x) = [x[1]^2 + x[2]^2 - 1]  # Example function: unit circle
fD(x) = [2x[1] 2x[2]]         # Jacobian matrix
x0 = [0.5, 0.5]               # Initial guess
```
Newton(f, fD, x0)
"""
function Newton(fun::Function, funD::Function, x0::Vector{Float64}; maxEval= 1000, εAccuracy= 1e-7)
    evalCount= 0
    improvementFound= true; direction= Vector{Float64}(undef, length(x0))
    xMin= x0; 
    fMin= fun(x0); 
    nfMin= fMin
    
    while (improvementFound || fMin> εAccuracy) && evalCount < maxEval   # run until no improvement found
        evalCount+= 1           # count evaluations of derivatives
        if improvementFound      # compute Newton's direction
            direction= vec(pinv(funD(xMin)) * fMin) 
            direction[isnan.(direction)] .= 0;
        else                    # nothing found: guess a direction
            direction= randn(length(x0))* (1e-7 + norm(direction))
        end
        α= 1.0; improvementFound= false
        while !improvementFound && (α > 0.6 || x0 ≠ xMin)   # never give up
            x0= xMin- α* direction  # handle NaN, Inf
            x0[isnan.(x0)] .= 0;
            if !all(isfinite.(x0))
              break
            end
           fx= fun(x0)         # function evaluation
           normf= fx
            if normf < nfMin    # improvement found
                xMin= x0; fMin= fx; nfMin= normf; improvementFound= true
            else                # half Newton step
                α/= 2           # ensure α will eventually be (exactly) 0
            end
        end
    end
    return (xMin= xMin, fMin= nfMin, evalCount= evalCount)
end

end #end of module 