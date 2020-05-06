# EnergyStatistics



*Energy Statistics (E-Statistics) for Julia.*

| **Build** | **Documentation**  |
|:---------:|:------------------:|
| [![Build Status][build-img]][build-url] | [![Docs][docs-dev-img]][docs-dev-url] |

This package allows to compute the distance covariance `dcov`,
distance variance `dvar` and distance correlation `dcor`.
See [this page][wiki] for more details and references.


#### Example

```julia
using EnergyStatistics

x = collect(-1:0.01:1)
y = map(x -> x^4 - x^2, x)

dcor(x, y) ≈ 0.374204050
```


[wiki]: https://en.wikipedia.org/wiki/Distance_correlation


[build-img]: https://travis-ci.org/pfarndt/EnergyStatistics.jl.svg?branch=master
[build-url]: https://travis-ci.org/pfarndt/EnergyStatistics.jl
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://pfarndt.github.io/EnergyStatistics.jl/dev
