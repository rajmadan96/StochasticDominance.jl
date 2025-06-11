```@meta
CurrentModule = StochasticDominance
```

# Introduction to Stochastic Dominance

Stochastic dominance plays a crucial role in decision-making under uncertainty and quantitative finance. It provides a powerful method for comparing random variables using their distribution functions. This concept establishes a structured framework to evaluate the superiority of one investment, policy, or strategy over another in uncertain environments.

By utilizing cumulative distribution functions, stochastic dominance allows decision-makers to assess choices without assuming specific utility functions. This generalization leads to a mathematically rigorous approach to optimization and risk management. The precision and depth of stochastic dominance make it an indispensable tool for analyzing complex probabilistic systems. We refer to [stochastic dominance wiki](https://en.wikipedia.org/wiki/Stochastic_dominance) for a quick overview. 


## Stochastic Dominance of Order $p$

Let $p \in [1,\infty)$ and let $X, Y \in L^p$ be random variables. We say that $X$ is dominated by $Y$ in the stochastic order of $p$, denoted as

```math
X \preccurlyeq^{(p)} Y
```

if the condition holds

```math
\mathbb{E} (t - X)_+^{p-1} \geq \mathbb{E} (t - Y)_+^{p-1}, \quad \text{for all } t \in \mathbb{R}.
```
This definition provides a precise way to compare random variables under higher order stochastic dominance criteria.
From a portfolio optimization perspective, the feasible random variable $Y$ in the above formulation dominates benchmark variables $X$. 

For other equivalent formulations and technical details, we refer to [Dentcheva2024](@citet).

## Optimal Portfolio Scaling under Stochastic Dominance of Order $p$

Given a **fixed benchmark** random vector $\xi_0 \in \mathbb{R}^{n}$ and a **portfolio return** random vector(s) ($d$ vectors) $\xi \in \mathbb{R}^{d \times n}$, we aim to determine an **optimal scaling factor** $x > 0$ ($n$ dimensional simplex vecto) such that the scaled portfolio $x^{\top}\xi$ dominates $\xi_0$ under stochastic dominance of order of $p$. Here, we have $d$ assets and $n$ scenarios.

### Goal:
Find $x > 0$ such that:

```math
\xi_0 \preccurlyeq^{(p)} x^{\top}\xi \quad \text{or} \quad \mathbb{E}(t - \xi_0)_+^{p-1} \geq \mathbb{E}(t - x^{\top}\xi )_+^{p-1}, \quad \forall t \in \mathbb{R}.
```

### Optimization Variants:

1. **Maximize Expected Return under Dominance Constraint**:
   ```math
   \max_{x > 0} \mathbb{E} x^{\top}\xi \quad \text{subject to} \quad \xi_0 \preccurlyeq^{(p)} x^{\top}\xi.
   ```

2. **Minimize Risk (e.g., Variance) under Dominance Constraint**:
   ```math
   \min_{x > 0} \text{Var}(x^{\top}\xi) \quad \text{subject to} \quad \xi_0 \preccurlyeq^{(p)} x^{\top}\xi.
   ```

### Interpretation:
- The condition $\xi_0 \preccurlyeq^{(p)} x^{\top}\xi$ ensures that for all thresholds $t$, the **$p^{\text{th}}$  partial moment** of $\xi_0$ is at least that of the scaled portfolio $x^{\top}\xi$.
- This allows the investor to **scale the portfolio** to either increase returns or reduce risk while maintaining or improving the performance relative to the benchmark under **higher-order stochastic dominance**.

### Practical Use:
This approach is especially useful in **portfolio optimization** where a benchmark (e.g., market index) is fixed, and the investor seeks to optimally scale their portfolio to satisfy risk-return tradeoffs under rigorous dominance conditions.

This package provides tools to analyze stochastic dominance, enabling users to apply these concepts effectively in decision-making and financial optimization. The mathematical background of this package was discussed in [Lakshmanan2025](@citet).

# References

```@bibliography
*
```



