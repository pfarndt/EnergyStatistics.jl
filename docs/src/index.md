# EnergyStatistics.jl


In statistics distance correlation or distance covariance is a measure of dependence between two paired random vectors. The population distance correlation coefficient is zero if and only if the random vectors are independent. Thus, distance correlation measures both linear and nonlinear association between two random vectors. This is in contrast to Pearson's correlation, which can only detect linear association between two random variables. See [here](https://en.wikipedia.org/wiki/Distance_correlation) for references and more details.


## Installation
This package can be installed using the Julia package manager. From the Julia REPL, type `]`
to enter the Pkg REPL mode and run

```
pkg> add EnergyStatistics
```

## General Usage
Given two random vectors `x` and `y` the distance correlation `dcor` can
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


Function to compute the distance covariance `dcov` and distance variance `dvar` are also supplied.



## Advanced Usage

The computation of the 'distance matrices' is computationally expansive.
In cases where you want to compute several distance correlation
keeping one of vectors fixed you may consider storing the distance matrices
for future use.

For example
```julia
A = dcenter!(DistanceMatrix(x))
B = dcenter!(DistanceMatrix(y))
C = dcenter!(DistanceMatrix(z))

dcor(A, B)
dcor(A, C)
```
will run  much faster than
```julia
dcor(x, y)
dcor(x, z)
```
since the distance matrix for the vector `x` has to be compute only once.

You can also construct distance matrices using other distance measures
than the (default) `abs`.

```julia
AA = dcenter!(DistanceMatrix(Float64, x, abs2))
```

In place double centering `dcenter!` and U-centering `ucenter!` functions are
available.


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
EnergyStatistics.dcor(A::EnergyStatistics.DistanceMatrix{T}, B::EnergyStatistics.DistanceMatrix{T},
        dvarA::T = dvar(A), dvarB::T = dvar(B) ) where {T <: Real}
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
