```@meta
CurrentModule = StochasticDominance
```

# StochasticDominance.jl

We present the `StochasticDominance.jl` package, which provides tools for analyzing higher order stochastic dominance—a method used to establish a partial order between random variables.

## Brief overview of stochastic dominance

Stochastic dominance is a concept used to compare decision alternatives based on their cumulative risk profiles, ensuring that one alternative is preferred over another without any trade-offs across values.
In the context of portfolio optimization, given a benchmark asset and a portfolio of assets, we seek an optimal allocation that maximizes a chosen objective (e.g., maximizing returns) while satisfying (higher-order) stochastic dominance constraints. 

## Main features of the package

The `StochasticDominance.jl` package provides tools to:

1. **Verify higher-order stochastic dominance** – Check whether a given portfolio satisfies higher-order dominance criteria relative to a benchmark asset.
2. **Determine the optimal allocation** – Find the asset allocation that maximizes a chosen objective (e.g., maximizing returns) while adhering to stochastic dominance constraints. It supports two key objective functions: maximizing expected returns and minimizing higher-order risk measures.  

## Installation

The package `StochasticDominance.jl` can be installed in Julia REPL as follows:

```julia
julia> using Pkg
julia> Pkg.add("StochasticDominance")
julia> using StochasticDominance
```

Once you have installed `StochasticDominance.jl`, we recommend going through the tutorials from beginning to end to understand how to use the package to verify stochastic dominance and determine the optimal allocation between a benchmark asset and a portfolio.

## Important functions in the package

The `StochasticDominance.jl` package provides several important functions, which are explained in detail in the tutorials section.
Here, we provide a brief overview of the functions.

1. `verify_dominance`: This function checks whether the given benchmark asset, represented as the random variable $X$, and the weighted portfolio asset, represented as the random variable $Y$, exhibit a dominance relationship for the specified stochastic order. This means that $Y$ consistently yields preferable outcomes over $X$ in the specified stochastic order. 

2. `optimize_max_return_SD`: This function determines the optimal asset allocation that maximizes expected returns for a given stochastic order (`SDorder`). Additionally, using `optimize_max_return_SD(; plots=true)`, users can generate a pie chart displaying the optimal allocation in percentages, along with the maximized expected returns and benchmark returns. The function also includes the option `optimize_max_return_SD(; verbose=true)`,  which allows users to imprint the convergence (or dominance) of the numerical method.

3. `optimize_min_riskreturn_SD`: This function determines the optimal asset allocation by minimizing higher-order risk measures for a given stochastic order (`SDorder`) while also indicating whether dominance is achieved. Additionally, using `optimize_min_riskreturn_SD(; plots=true)`, users can generate a pie chart that visualizes the optimal allocation in percentages, along with the minimizing higher-order risk measure returns. The function also provides the option `optimize_min_riskreturn_SD(; verbose=true)`, allowing users to assess the convergence (or dominance) of the numerical algorithm.

