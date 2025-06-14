```@meta
CurrentModule = StochasticDominance
```

# Higher Order Stochastic Dominance Verification

This tutorial demonstrates the function `verify_dominance`, which validates the dominance relationship between $Y$ and $X$.

## Example 1
To verify if $Y$ stochastically dominates $X$ under the specified order `SDorder = 2`, we use the function `verify_dominance`

### Defining Random Variables $X$ and $Y$

Consider two discrete random variables, $X$ and $Y$. Their values and associated probabilities are:

- $Y = [3, 5, 7]$ with probabilities $p_Y = [0.3, 0.4, 0.3]$.
- $X = [2, 4, 6]$ with probabilities $p_X = [0.2, 0.5, 0.3]$.

These represent two probability distributions that we compare using stochastic dominance.

```julia
Y = [3, 5, 7]
X = [2, 4, 6]
p_Y = [0.3, 0.4, 0.3]
p_X = [0.2, 0.5, 0.3]
SDorder = 2
verify_dominance(Y, X,SDorder;p_Y, p_X)

```
**Output:**
```julia
Y dominates X in stochastic order 2
true
```

This function checks whether $Y$ stochastically dominates $X$ of order `SDorder`. The `verbose=true` option ensures detailed output, providing insights into the dominance verification process.


## Example 2

Consider two discrete random variables, $X$ and $Y$, with swapped roles. Their values and associated probabilities are:

- $X = [3, 5, 7]$ with probabilities $p_X = [0.3, 0.4, 0.3]$.
- $Y = [2, 4, 6]$ with probabilities $p_Y = [0.2, 0.5, 0.3]$.

These represent two probability distributions that we compare using stochastic dominance.

```julia
X = [3, 5, 7]
Y = [2, 4, 6]
p_X = [0.3, 0.4, 0.3]
p_Y = [0.2, 0.5, 0.3]
SDorder = 2
verify_dominance(Y,X,SDorder;p_Y,p_X,verbose=true)
```
**Output:**
```julia
Y doesn't dominates X in stochastic order 2
Y doesn't dominates X in stochastic order 2 with residual 1.6970562748477143
false
```

The verbose output indicates that with a residual value of 1.70, we can confirm that Y does not dominate X in the second-order stochastic dominance. However, in general, the user has control over the choice of residual (ε), and it can be selected as needed.

