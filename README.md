[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://rajmadan96.github.io/StochasticDominance.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://rajmadan96.github.io/StochasticDominance.jl/dev/)
[![Aqua QA](https://img.shields.io/badge/Aqua.jl-%F0%9F%8C%A2-aqua.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Coverage](https://codecov.io/gh/rajmadan96/StochasticDominance.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/rajmadan96/StochasticDominance.jl)

# StochasticDominance.jl

`Stochasticominance.jl` is a Julia package that offers tools for analyzing higher-order stochastic dominance, a technique used to define a partial order between random variables.

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

1. `VerifyDominance`: This function checks whether the given benchmark asset and the weighted portfolio assets exhibit a dominance relationship for the specified stochastic order. Note that the benchmark asset scenario and the weighted portfolio assets must have the same length for a valid comparison.

2. `StochasticDominanceMeanReturn`: This function determines the optimal asset allocation that maximizes expected returns for a given stochastic order (`SDorder`).
Additionally, using `StochasticDominanceMeanReturn(; plots=true)`, users can generate a pie chart displaying the optimal allocation in percentages, along with the maximized expected returns and benchmark returns. The function also includes the option `StochasticDominanceMeanReturn(; verbose=true)`, which allows users to evaluate the convergence (or dominance) quality.

3. `StochasticDominanceRiskMeasure`: This function determines the optimal asset allocation by minimizing higher-order risk measures for a given stochastic order (`SDorder`) while also indicating whether dominance is achieved. Additionally, using `StochasticDominanceRiskMeasure(; plots=true)`, users can generate a pie chart that visualizes the optimal allocation in percentages, along with the maximized expected returns and benchmark returns. The function also provides the option `StochasticDominanceRiskMeasure(; verbose=true)`, allowing users to assess the convergence (or dominance) quality.

