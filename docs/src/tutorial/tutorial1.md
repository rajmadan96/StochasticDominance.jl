```@meta
CurrentModule = StochasticDominance
```

# Introduction to Stochastic Dominance

Stochastic dominance plays a crucial role in decision-making under uncertainty and quantitative finance. It provides a powerful method for comparing random variables using their distribution functions. This concept establishes a structured framework to evaluate the superiority of one investment, policy, or strategy over another in uncertain environments.

By utilizing cumulative distribution functions, stochastic dominance allows decision-makers to assess choices without assuming specific utility functions. This generalization leads to a mathematically rigorous approach to optimization and risk management. The precision and depth of stochastic dominance make it an indispensable tool for analyzing complex probabilistic systems.

## Mathematical Framework

We formulate the stochastic dominance problem as follows
```math
\begin{aligned}
    \text{maximize } & \mathbb{E}\,Y \\
    \text{subject to } & X \preccurlyeq Y \text{ for every } X \in \mathcal{X}.
\end{aligned}
```

Here, the notation “$\preccurlyeq$” represents a partial order between random variables $X$ and $Y$, where $Y$ belongs to the set $\mathcal{X}$, which consists of benchmark random variables. More specifically, we consider general stochastic dominance relations of order $p \geq 1$, denoted by $\preccurlyeq^{(p)}$.

From a portfolio optimization perspective, the feasible random variable $Y$ in the above formulation dominates all benchmark variables $X$. Often, the set $\mathcal{X}$ contains a single random variable $X$, reducing the problem to finding a random variable $Y$ such that $X \preccurlyeq Y$ while optimizing the objective.

## Definition: Stochastic Dominance of Order $p$

Let $p \in [1,\infty)$ and let $X, Y \in L^p$ be random variables. We say that $X$ is dominated by $Y$ in the $p$th stochastic order, denoted as

```math
X \preccurlyeq^{(p)} Y
```

if the condition holds

```math
\mathbb{E} (t - X)_+^{p-1} \geq \mathbb{E} (t - Y)_+^{p-1}, \quad \text{for all } t \in \mathbb{R}.
```

Additionally, we define

```math
X \preccurlyeq^{(\infty)} Y \quad \text{if} \quad \operatorname{essinf} X \leq \operatorname{essinf} Y.
```

This definition provides a precise way to compare random variables under higher order stochastic dominance criteria.

For other equivalent formulations and technical details, we refer to [Dentcheva and Ruszczyński (2024)](https://doi.org/10.1007/s10589-015-9758-0).

The stochastic dominance of order $p$ and stochastic dominance with respect to the norm $\|\cdot\|_p$ are related by

```math
\preccurlyeq^{(p+1)} ~~~ \iff ~~~ \preccurlyeq^{\|\cdot\|_p}.
```
The reason for the apparent, unfortunate mismatch between $p$ and $p+1$ is historic. 
In other words, whenever we refer to stochastic order (`SDorder`), it means

```math
SDorder = p - 1.
```
 
This package provides tools to analyze stochastic dominance, enabling users to apply these concepts effectively in decision-making and financial optimization.

