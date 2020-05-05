# EnergyStatistics.jl


In statistics distance correlation or distance covariance is a measure of dependence between two paired random vectors. The population distance correlation coefficient is zero if and only if the random vectors are independent. Thus, distance correlation measures both linear and nonlinear association between two random vectors. This is in contrast to Pearson's correlation, which can only detect linear association between two random variables. See [here](https://en.wikipedia.org/wiki/Distance_correlation) for references and more details.


## Installation
This package can be installed using the Julia package manager. From the Julia REPL, type `]`
to enter the Pkg REPL mode and run

```
pkg> add EnergyStatistics
```

## General Usage
Given two vectors `x` and `y` the distance correlation `dcor` can
simply computed:

```julia
using EnergyStatistics

x = collect(-1:0.01:1)
y = map(x -> x^4 - x^2, x)

dcor(x, y) ≈ 0.374204050
```

These two vectors are clearly associated. However, their (Pearson)
correlation coefficient vanishes suggesting that they are independent.
The finite distance correlation `dcor`
reveals their __non-linear__ association.


Function to compute the distance covariance `dcov` and distance variance `dvar`
are also supplied.



## Advanced Usage

The computation of a 'DistanceMatrix' is computationally expansive.
Especially the computation of `n(n-1)/2` pairwise distances
for vectors of length `n` and the subsequent
centering of the distance matrix take time and memory.
In cases where one wants to compute several distance correlations and
keep intermediate results of the distance computations and centering
one can do so. For example:

```julia
Dx = EnergyStatistics.dcenter!(EnergyStatistics.DistanceMatrix(x))
Dy = EnergyStatistics.dcenter!(EnergyStatistics.DistanceMatrix(y))
Dz = EnergyStatistics.dcenter!(EnergyStatistics.DistanceMatrix(z))

dcor_xy = dcor(Dx, Dy)
dcor_xz = dcor(Dx, Dz)
```
will run faster than

```julia
dcor_xy = dcor(x, y)
dcor_xz = dcor(x, z)
```
since the distance matrix `Dx` for the vector `x` is only computed once.


You can also construct distance matrices using other distance measures
than the (default) `abs`.

```julia
AA = EnergyStatistics.dcenter!(EnergyStatistics.DistanceMatrix(Float64, x, abs2))
```


Instead of double centering via `dcenter!` one may also
use U-centering via the `ucenter!` function.


## References

See the [wikipedia page](https://en.wikipedia.org/wiki/Distance_correlation) for references and more details.



## Functions


```@docs
EnergyStatistics.dcor(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Real}
```

```@docs
EnergyStatistics.dcov(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Real}
```

```@docs
EnergyStatistics.dvar(x::AbstractVector{T}) where {T <: Real}
```




```@docs
EnergyStatistics.DistanceMatrix(x::AbstractVector{T}, dist = abs ) where {T}
```

```@docs
EnergyStatistics.dcenter!(A::EnergyStatistics.DistanceMatrix{T}) where {T <: Real}
```

```@docs
EnergyStatistics.ucenter!(A::EnergyStatistics.DistanceMatrix{T}) where {T <: Real}
```

```@docs
EnergyStatistics.dcor(A::EnergyStatistics.DistanceMatrix{T}, B::EnergyStatistics.DistanceMatrix{T}) where {T <: Real}
```

```@docs
EnergyStatistics.dcov(A::EnergyStatistics.DistanceMatrix{T}, B::EnergyStatistics.DistanceMatrix{T}) where {T <: Real}
```

```@docs
EnergyStatistics.dvar(A::EnergyStatistics.DistanceMatrix{T}) where {T <: Real}
```


## Index

```@index
```
