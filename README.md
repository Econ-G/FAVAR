## Acknowledgments

Huge thanks to **Yuchen** for her invaluable support throughout this project. Your insights, patience, and timely feedback made all the difference in bringing this work to life!

## Data
We use 186 Canadian
macroeconomic monthly time series and the monthly average federal funds rate from January 1997 to November 2024.

## Model

This paper uses the FAVAR approach, developed by (Bernanke et al., 2005)
to study the spillover effect of US monetary policy on the Canadian
economy. Let $Y_t$ be a M $\times$ 1 vector of observable economic
series and F be a K $\times$ 1 vector of unobservable factors, where M
is the number of series in $Y_t$ and K the number of factors. In several
cases, the dynamics of the economy are not fully represented by $Y_t$.
Factors $F_t$ are designed to capture additional information. The joint dynamics of (F_t, Y_t) are assumed to follow this transition equation:


```math

  \begin{bmatrix} 
    F_t \\ 
    Y_t 
  \end{bmatrix}
  = \Phi(L)\,
    \begin{bmatrix} 
      F_{t-1} \\ 
      Y_{t-1} 
    \end{bmatrix}
    + v_t
  \space\space\space\space\space(1)
```


where $\Phi$ is a lag polynomial of finite order $p$ (the
number of lags), and $v_t$ is the error term with mean zero and
covariance matrix $Q$. The equation (1)
cannot be directly estimated because $F_t$ is unobservable. The factor
can be interpreted as representing the underlying forces that can
potentially affect many economic variables. We can then uncover those
factors reservedly. Here is a simplified example:

```math

F_{t}^{\mathrm{Business\ Cycle}}
  = \alpha_1\,\mathrm{GDP}
  + \alpha_2\,\mathrm{Unemployment}
  + \alpha_3\,\mathrm{Interest\ Rate}
  + \alpha_4\,\mathrm{Consumption}
  + \varepsilon_t

```

where the $F_t$ is an unobservable variable that represents the business
cycle. We include explanatory variables, as the observable variables
representing macroeconomics that move along the business cycle. If we
assume that this regression is a valid relationship, we can study the
common variance shared by independent variables to estimate the
unobservable latent factors. We set $X_t$ as an N $\times$ 1 matrix,
which contains an extensive amount of data sets. The $X_t$ is correlated
with $F_t$ and $Y_t$ by the following equation:

```math

X_t = \Lambda^f F_t + \Lambda^y Y_t + u_t
\space\space\space\space\space(2)
```

where $\Lambda^f$ and $\Lambda^y$ are the factor loadings relating $F_t$
and $Y_t$ to the data in $X_t$. The $\Lambda^f$ is an N $\times$ M
matrix and $\Lambda^y$ is an N $\times$ M matrix. $u_t$ is the vector of
error terms, which are assumed to have zero mean and to be uncorrelated
with the elements of $F_t$ and $Y_t$. Then we substitute $F_t$ with
$\hat{F}$ in two steps. First, we extract principal components $\hat{C}$
from $X_t$ (186 Canadian macro series) and $Y_t$ Federal Funds Rate. All
principal components are normalized to have unit variances. Second, to
ensure the identification of the VAR model, the latent factors cannot
respond contemporaneously to monetary policy innovations. Hence, to
implement this identification scheme, we separate all Canadian variables
into two categories: "fast-moving" and "slow-moving" variables. A
variable is classified as "slow-moving" if it shows no contemporaneous
response to monetary shocks, for example, the Consumer Price Index. The
"fast-moving" variables are very sensitive to monetary shocks, such as
stock prices and exchange rates. The classification of 186 macro
variables is provided in the data section. The "slow-moving" factors
$\hat{F^s_t}$ are estimated as the principal components of the
"slow-moving" variables. And they do not reacting contemporaneously to
the monetary shock. We can then regress the regression 
to obtain $b_{F^s}$ and $b_Y$. The latent
factors $\hat{F_t}$ then can be constructed from

```math
\hat{C}_t = b_{F^s} \hat{F}_t^s + b_y Y_t + e_t
\space\space\space\space\space(3)
```


