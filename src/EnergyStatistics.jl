module EnergyStatistics

export dcov, dvar, dcor


using LinearAlgebra

const DistanceMatrix{T} = Matrix{T}


function DistanceMatrix(::Type{F}, x::AbstractVector{T}, dist) where {F <: Real, T}
    n = length(x)
    A = DistanceMatrix{F}(undef, n, n)

    d0 = zero(F)
    for i in 1:n
        @inbounds xi = x[i]
        @inbounds A[i, i] = d0
        for j in (i+1):n
            @inbounds dij = dist(xi - x[j])
            @inbounds A[i, j] = dij
            @inbounds A[j, i] = dij
        end
    end
    A
end




"""
    DistanceMatrix(x::AbstractVector{T}, dist = abs) where {T}

Computes the matrix of pairwise distance of `x`. The distance measure `dist`
is `abs` as default.

```jldoctest
using EnergyStatistics
x = [1.0, 2.0]
EnergyStatistics.DistanceMatrix(x)

# output

2×2 Array{Float64,2}:
 0.0  1.0
 1.0  0.0
```
"""
DistanceMatrix(x::AbstractVector{T}, dist = abs) where {T} = DistanceMatrix(Float64, x, dist)

"""
    dcenter!(A::DistanceMatrix{T}) where {T <: Real}

Computes the double centered matrix of `A` in place.
"""
function dcenter!(A::DistanceMatrix{T}) where {T <: Real}
    n = LinearAlgebra.checksquare(A)

    # double center of symmetric matrix
    s = sum(A, dims = 1) ./ n
    ss = sum(s) / n

    for i in 1:n
        @inbounds A[i, i] = A[i, i] - 2 * s[i] + ss
        for j in (i+1):n
            @inbounds aij = A[i, j] - s[i] - s[j] + ss
            @inbounds A[i, j] = aij
            @inbounds A[j, i] = aij
        end
    end
    A
end


"""
    ucenter!(A::DistanceMatrix{T}) where {T <: Real}

Computes the u-centered matrix of `A` in place.
"""
function ucenter!(A::DistanceMatrix{T}) where {T <: Real}
    n = LinearAlgebra.checksquare(A)
    n > 2 || throw("Matrix too small for u-centering")

    # double center of symmetric matrix
    s = sum(A, dims = 1)
    ss = sum(s)

    s =  s ./ (n-2)
    ss = ss / (n-1)*(n-2)

    for i in 1:n
        @inbounds A[i, i] = zero(T)
        for j in (i+1):n
            @inbounds aij = A[i, j] - s[i] - s[j] + ss
            @inbounds A[i, j] = aij
            @inbounds A[j, i] = aij
        end
    end
    A
end


"""
    dcov(A::DistanceMatrix{T}, B::DistanceMatrix{T}) where {T <: Real}

Computes the distance covariance of two doubly-centered distance matrices.
"""
function dcov(A::DistanceMatrix{T}, B::DistanceMatrix{T}) where {T <: Real}
    ns = LinearAlgebra.checksquare(A, B)
    ns[1] == ns[2] || throw(DimensionMismatch("matrices must be of same size"))
    n = ns[1]

    s = zero(T)
    for i in 1:n
        @inbounds s += A[i, i] * B[i, i]
        for j in i+1:n
            @inbounds s += 2 * A[j, i] * B[j, i]
        end
    end
    sqrt(s / n^2)
end


"""
    dvar(A::DistanceMatrix{T}) where {T <: Real}

Computes the distance variance of a doubly-centered distance matrix.
"""
function dvar(A::DistanceMatrix{T}) where {T <: Real}
    n = LinearAlgebra.checksquare(A)

    s = zero(T)
    for i in 1:n
        @inbounds s += A[i, i]^2
        for j in i+1:n
            @inbounds s += 2 * A[j, i]^2
        end
    end
    sqrt(s / n^2)
end


"""
    dcor(A::DistanceMatrix{T}, B::DistanceMatrix{T}, dvarA::T = dvar(A), dvarB::T = dvar(B))

Computes the distance correlation of two DistanceMatrices `A` and `B`.
"""
function dcor(A::DistanceMatrix{T}, B::DistanceMatrix{T},
        dvarA::T = dvar(A), dvarB::T = dvar(B) ) where {T <: Real}
    dvarAdvarb = dvarA * dvarB
    iszero(dvarAdvarb) && return zero(T)
    dcov(A, B) / sqrt(dvarAdvarb)
end



"""
    dcov(x::AbstractVector{T}, y::AbstractVector{T}) where T <: Real

Computes the distance covariance of samples `x` and `y`.
"""
function dcov(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Real}
    length(x) == length(y) || throw(DimensionMismatch("vectors must be of same length"))
    A = dcenter!(DistanceMatrix(x))
    B = dcenter!(DistanceMatrix(y))
    dcov(A, B)
end


"""
    dvar(x::AbstractVector{T}) where T <: Real

Computes the distance variance of a sample `x`.
"""
function dvar(x::AbstractVector{T}) where {T <: Real}
    A = dcenter!(DistanceMatrix(x))
    dvar(A)
end


"""
    dcor(x::AbstractVector{T}, y::AbstractVector{T}) where T <: Real

Computes the distance correlation of samples `x` and `y`.

```jldoctest
using EnergyStatistics
x = collect(-1:0.01:1)
y = @. x^4 - x^2
dcor(x, y)

# output

0.3742040504583154
```
"""
function dcor(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Real}
    length(x) == length(y) || throw(DimensionMismatch("vectors must be of same length"))
    A = dcenter!(DistanceMatrix(x))
    B = dcenter!(DistanceMatrix(y))
    dcor(A, B)
end

end # module
