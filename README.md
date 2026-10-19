# Persistance of productivity shocks in stochastic growth models

A MATLAB program that shows the effect of changing the persistence of productivity shocks in a simple neoclassical growth model where these shocks come in the form of a finite-state Markov process. Evidently, increasing the persistence of productivity states smoothes consumption and investment, while decreasing it (i.e. increasing volatility) will result in more volatile consumption. 

As the model also shows, higher shock persistence broadens the capital distribution (e.g., it is easier for the agent to find himself with very high/low capital levels). On the contrary, more volatility in the productivity states concentrates the distribution at the mean. The borrowing constraint becomes more relevant as persistence increases: since the agent is more likely to spend a large amount of periods in a low-productivity state, the fraction of periods spent at the constraint is higher and so is the probability of hitting the constraint.

## The Model

The planner solves

$$
\max_{\{c_t\}_{t=0}^{\infty}} \; \mathbb{E}_0 \sum_{t=0}^{\infty} \beta^t \ln(c_t)
$$

subject to the resource constraint

$$
c_t + k_{t+1} = A_t k_t^{\alpha} + (1-\delta)k_t
$$

and the borrowing constraint

$$
k_{t+1} \geq \underline{k}
$$

Productivity $A_t$ takes two values, $A_H$ and $A_L$, and follows a first-order Markov process with transition matrix $\Pi$.

### Parameters

| Parameter | Description | Value |
|---|---|---|
| $\beta$ | Discount factor | 0.95 |
| $\alpha$ | Capital share of income | 0.33 |
| $\delta$ | Depreciation rate | 0.04 |
| $A_H$, $A_L$ | High / low TFP | 1.1 / 0.9 |
| $\underline{k}$ | Minimum capital level | 6.2 |

### Transition Matrices

Default Scenario | More Persistence | More Volatility |
|---|---|---|
| $$\begin{bmatrix} 0.51 & 0.49 \\ 0.49 & 0.51 \end{bmatrix}$$ | 4$\begin{bmatrix} 0.60 & 0.40 \\ 0.40 & 0.60 \end{bmatrix}$$ | $$\begin{bmatrix} 0.40 & 0.60 \\ 0.60 & 0.40 \end{bmatrix}$$ |

## Output

The script reports:
- Mean and standard deviation of capital in the stationary distribution
- Probability that the borrowing constraint binds
- Fraction of simulated periods spent at the constraint

and produces five figures:
1. Value functions for each TFP state
2. Policy functions with the 45° line and the constrained region shaded
3. Stationary distribution of capital
4. Simulated consumption, investment and productivity shocks
5. Simulated capital path, with periods at the constraint highlighted

To make comparison between scenarios easier, each figure can be toggled on/off at the very top of the script. This ever so slightly increases performance and prevents 15 figures from spamming your screen.

## Usage

Open MATLAB and run:

```matlab
persistence_volatility_example
```

That's it, no gimmicks.

## Author

Tomás Roda Faria Rodrigues
