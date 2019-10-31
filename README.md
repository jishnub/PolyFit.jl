# PolyFit

<!---
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jishnub.github.io/PolyFit.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jishnub.github.io/PolyFit.jl/dev)
-->
[![Build Status](https://travis-ci.com/jishnub/PolyFit.jl.svg?branch=master)](https://travis-ci.com/jishnub/PolyFit.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jishnub/PolyFit.jl?svg=true)](https://ci.appveyor.com/project/jishnub/PolyFit-jl)
[![Codecov](https://codecov.io/gh/jishnub/PolyFit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jishnub/PolyFit.jl)
[![Coveralls](https://coveralls.io/repos/github/jishnub/PolyFit.jl/badge.svg?branch=master)](https://coveralls.io/github/jishnub/PolyFit.jl?branch=master)
[![Build Status](https://api.cirrus-ci.com/github/jishnub/PolyFit.jl.svg)](https://cirrus-ci.com/github/jishnub/PolyFit.jl)

This package extends [Polynomials.jl](https://github.com/JuliaMath/Polynomials.jl) and provides methods to fit a polynomial to data with error bars. The solution is obtained through an L2 minimization, where the covariance matrix is used to define the misfit that is minimized.

We solve the system 

![equation](https://latex.codecogs.com/gif.latex?%5Cleft%5C%7B%20%5Cbegin%7Barray%7D%7Bc%7D%20y_%7B1%7D%5C%5C%20y_%7B2%7D%5C%5C%20%5Cvdots%5C%5C%20y_%7Bm%7D%20%5Cend%7Barray%7D%5Cright%5C%7D%20%3D%5Cleft%5B%5Cbegin%7Barray%7D%7Bccccc%7D%201%20%26%20x_%7B1%7D%20%26%20x_%7B1%7D%5E%7B2%7D%20%26%20%5Ccdots%20%26%20x_%7B1%7D%5E%7Bn-1%7D%5C%5C%201%20%26%20x_%7B2%7D%20%26%20x_%7B2%7D%5E%7B2%7D%20%26%20%5Ccdots%20%26%20x_%7B2%7D%5E%7Bn-1%7D%5C%5C%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cddots%20%26%20%5Cvdots%5C%5C%201%20%26%20x_%7Bm%7D%20%26%20x_%7Bm%7D%5E%7B2%7D%20%26%20%5Ccdots%20%26%20x_%7Bm%7D%5E%7Bn-1%7D%20%5Cend%7Barray%7D%5Cright%5D%5Cleft%5C%7B%20%5Cbegin%7Barray%7D%7Bc%7D%20a_%7B0%7D%5C%5C%20a_%7B1%7D%5C%5C%20%5Ccdots%5C%5C%20a_%7Bn-1%7D%20%5Cend%7Barray%7D%5Cright%5C%7D%20)

We can express this as 

![equation](https://latex.codecogs.com/gif.latex?%5Cmathbf%7By%7D%3DX%5Cmathbf%7Ba%7D)

Solving the equation system is equivalent to minimizing the misfit

![equation](https://latex.codecogs.com/gif.latex?%5Cchi%5E%7B2%7D%3D%5Cleft%5C%7B%20X%5Cmathbf%7Ba%7D-%5Cmathbf%7By%7D%5Cright%5C%7D%20%5E%7BT%7DC%5E%7B-1%7D%5Cleft%5C%7B%20X%5Cmathbf%7Ba%7D-%5Cmathbf%7By%7D%5Cright%5C%7D,)

where `C` is the covariance matrix. The diagonal elements of the covariance matrix are the variances in `y`. Fitting data without error bars is equivalent to setting `C = σ^-2 I` for a constant `σ`.

# Installation

Install the package using 

```julia
pkg> add https://github.com/jishnub/PolyFit.jl

julia> using PolyFit
```

# Usage

This package adds methods to `polyfit`, and the usage is similar to `Polynomials.jl`

```julia
julia> x=1:50;y=x.^2;

julia> polyfit(x,y,2)
Poly(3.1065291710470317e-14 - 2.6679598459471846e-15*x + 1.0000000000000002*x^2)
```

The extra functionality this adds is the ability to pass error bars on `y` to the fit. This can be done either by passing error bars `σ(y)` or a covariance matrix `C(yi,yj)`.

## Error bars

```julia
julia> σ=10rand(length(y));

julia> polyfit(x,y,2,σ)
Poly(2.643961948331425e-14 - 2.7668654490519844e-15*x + 1.0000000000000002*x^2)
```

## Covariance Matrix

A diagonal covariance matrix is equiavlent to passing error bars. The diagonal elements are the variances.

```julia
julia> polyfit(x,y,2,Diagonal(σ.^2))
Poly(2.6426223525629683e-14 - 2.7636271626091054e-15*x + 1.0000000000000002*x^2)
```

A non-diagonal covariance matrix has to be positive-definite

```julia
julia> C = Symmetric(Diagonal(σ.^2) + 0.1*rand(length(y),length(y)));

julia> isposdef(C)
true

julia> polyfit(x,y,2,C)
Poly(5.517934888038756e-14 + 5.126632798188113e-15*x + 0.9999999999999998*x^2)
```

## SpecialMatrices

This package uses [SpecialMatrices.jl](https://github.com/JuliaMatrices/SpecialMatrices.jl) to solve Van der Monde matrices, so the special case of fitting a degree-n polynomial to n+1 points without error bars is faster than `Polynomials.jl`. This method can be called using the special tag `SquareSystem()`.

```julia
julia> @btime polyfit($x,$y,length($x)-1);
  168.393 μs (59 allocations: 65.02 KiB)

julia> @btime polyfit(SquareSystem(),$x,$y,length($x)-1);
  2.071 μs (4 allocations: 1.11 KiB)

julia> polyfit(x,y,length(x)-1)≈polyfit(SquareSystem(),x,y,length(x)-1)
true
```