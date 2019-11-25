# Distance Correlation



*Distance Correlation for Julia.*

| **Documentation**  |
|:------------------:|
| [![][docs-dev-img]][docs-dev-url] |

This package allows to compute the distance covariance `dcov`, distance variance `dvar` and distance correlation `dcor`.
See [this page][wiki] for more details and references.


#### Example

```julia
using EnergyStatistics

x = collect(-1:0.01:1)
y = map(x -> x^4 - x^2, x)

dcor(x, y) ≈ 0.374204050
```


[wiki]: https://en.wikipedia.org/wiki/Distance_correlation

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://pfarndt.github.io/EnergyStatistics.jl/dev
