module PSOOpt

using LinearAlgebra, Random
Random.seed!(42)

export PSO

"""
    PSO(fun::Function, x0::Vector{Float64}; ω = 0.5, φp = 1.5, φg = 1.5, 
        n_particles = 30, maxEval = 10000, threshold = 1e-9)

Performs Particle Swarm Optimization (PSO) to minimize a given objective function.

# Arguments
- `fun::Function`: The objective function to be minimized. It should accept a vector of variables as input and return a scalar or vector.
- `x0::Vector{Float64}`: An initial guess for the solution. The dimensionality of the optimization is inferred from the length of this vector.

# Keyword Arguments
- `ω::Float64 = 0.5`: Inertia weight controlling the contribution of the particle's previous velocity.
- `φp::Float64 = 1.5`: Cognitive coefficient that scales the particle's tendency to move toward its personal best.
- `φg::Float64 = 1.5`: Social coefficient that scales the particle's tendency to move toward the global best.
- `n_particles::Int = 30`: Number of particles in the swarm.
- `maxEval::Int = 10000`: Maximum number of iterations for the optimization process.
- `threshold::Float64 = 1e-9`: Convergence threshold; the optimization stops when the global best score falls below this value.

# Returns
A tuple containing:
- `xMin::Vector{Float64}`: The position of the global best solution found by the swarm.
- `fMin::Float64`: The value of the objective function at `xMin`.
```
"""

# Define the PSO function
function PSO(fun::Function,x0::Vector{Float64};ω = .5,φp = 1.5,φg = 1.5,n_particles = 30,maxEval = 1000)
    # PSO parameters
        
    dim = length(x0) # Number of variables
    # Initialize particle positions, velocities, and personal/global bests
    positions = [rand(dim) for _ in 1:n_particles]
    velocities = [rand(dim) for _ in 1:n_particles] 
    personal_best_positions = copy(positions)
    personal_best_scores = [fun(q) for q in positions]
    global_best_position = personal_best_positions[argmin(personal_best_scores)]
    global_best_score = minimum(personal_best_scores)

    # Particle Swarm Optimization loop
    for iter in 1:maxEval
        for i in 1:n_particles
            # Update velocity
            velocities[i] = ω * velocities[i] +
                            φp * rand(dim) .* (personal_best_positions[i] - positions[i]) +
                            φg * rand(dim) .* (global_best_position - positions[i])

            # Update position
            positions[i] += velocities[i]
            
            # Evaluate fitness
            current_score = (fun(positions[i]))
            if current_score < personal_best_scores[i]
                personal_best_positions[i] = positions[i]
                personal_best_scores[i] = current_score
            end

            # Update global best
            if current_score < global_best_score
                global_best_position = positions[i]
                global_best_score = current_score
            end
        end

        if iter === maxEval
           break
        end
    end

    return (xMin= global_best_position, fMin= global_best_score)
end


end # end of module