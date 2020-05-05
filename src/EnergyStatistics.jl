module EnergyStatistics

export dcov, dvar, dcor



# Distance Marices are symmetric and only the upper triagonal needs to be
# stored. The memory layout is optimized for computation od dcov/dvar/dcor.
mutable struct DistanceMatrix{T} <: AbstractMatrix{T}
    n::Int
    data::Vector{T}
    centering::Symbol
    dvar::Union{T, Nothing}
end

DistanceMatrix{T}(n::Int) where T = DistanceMatrix{T}(n, Vector{T}(undef, n * (n+1) ÷ 2), :none, nothing)

function dataindex(n::Int, i::Int, k::Int)
    Δ = abs(i - k)
    Δ * n - (Δ * (Δ - 1)) ÷ 2 + min(i, k)
end

Base.getindex(m::DistanceMatrix{T},i, k) where T = getindex(m.data, dataindex(m.n, i, k))
Base.setindex!(m::DistanceMatrix{T}, X, i, k) where T = setindex!(m.data, X, dataindex(m.n, i, k))
Base.size(m::DistanceMatrix{T}) where T = (m.n, m.n)
Base.eltype(m::DistanceMatrix{T}) where T = T


# To access elements some Iterator along data vector, i.e. diagonals
mutable struct DistanceMatrixDataIter{T}
    A::DistanceMatrix{T}
    i::Int
    j::Int
end

DistanceMatrixDataIter{T}(A::DistanceMatrix{T}) where T = DistanceMatrixDataIter{T}(A, 0, 0)
DistanceMatrixDataIter(A::DistanceMatrix{T}) where T = DistanceMatrixDataIter{T}(A)
Base.eltype(::Type{DistanceMatrixDataIter{T}}) where T = Tuple{Int, Int, Int ,T}
Base.length(dmi::DistanceMatrixDataIter{T}) where T = length(dmi.A.data)

function Base.iterate(dmi::DistanceMatrixDataIter{T}, state = 1) where T
    dmi.i += 1
    dmi.j += 1
    if dmi.j > dmi.A.n
        dmi.i == 2 && return(nothing)
        dmi.j = dmi.j - dmi.i + 2
        dmi.i = 1
    end
    ((dmi.i, dmi.j, state, dmi.A.data[state]), state+1)
end

# Iterator along rows
mutable struct DistanceMatrixRowIter{T}
    A::DistanceMatrix{T}
    i::Int
    j::Int
end

DistanceMatrixRowIter{T}(A::DistanceMatrix{T}) where T = DistanceMatrixRowIter{T}(A, 1, 0)
DistanceMatrixRowIter(A::DistanceMatrix{T}) where T = DistanceMatrixRowIter{T}(A)
Base.eltype(::Type{DistanceMatrixRowIter{T}}) where T = Tuple{Int, Int, Int ,T}
Base.length(dmi::DistanceMatrixRowIter{T}) where T = length(dmi.A.data)

function Base.iterate(dmi::DistanceMatrixRowIter{T}, state = 1) where T
    dmi.j += 1
    if dmi.j > dmi.A.n
        dmi.i ==  dmi.A.n && return(nothing)
        dmi.i += 1
        dmi.j = dmi.i
    end
    di = dataindex(dmi.A.n, dmi.i, dmi.j)
    ((dmi.i, dmi.j, di, dmi.A.data[di]), state+1)
end




function DistanceMatrix(::Type{F}, x::AbstractVector{T}, dist) where {F <: Real, T}
    n = length(x)
    A = DistanceMatrix{F}(n)

    for (i, j, di, v) in DistanceMatrixRowIter(A)
        if i == j
            @inbounds A.data[di] = zero(F)
        else
            @inbounds A.data[di] = dist(x[i] - x[j])
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

2×2 EnergyStatistics.DistanceMatrix{Float64}:
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
    n = size(A, 1)

    # double center of symmetric matrix
    si = zeros(T, n)
    sj = zeros(T, n)
    for (i, j, di, v) in DistanceMatrixDataIter(A)
        if i == j
            si[i] += v
        else
            si[i] += v
            sj[j] += v
        end
    end
    s = (si .+ sj) ./ n
    ss = sum(s) / n

    for (i, j, di, v) in DistanceMatrixDataIter(A)
        A.data[di] = A.data[di] - s[i] - s[j] + ss
    end
    A.centering = :dcenter
    A.dvar = nothing
    A
end


"""
    ucenter!(A::DistanceMatrix{T}) where {T <: Real}

Computes the u-centered matrix of `A` in place.
"""
function ucenter!(A::DistanceMatrix{T}) where {T <: Real}
    n = size(A, 1)
    n > 2 || throw("Matrix too small for u-centering")


    # double center of symmetric matrix
    s = zeros(T, n)
    sj = zeros(T, n)
    for (i, j, di, v) in DistanceMatrixDataIter(A)
        if i == j
            s[i] += v
        else
            s[i] += v
            sj[j] += v
        end
    end
    s = s .+ sj
    ss = sum(s)
    s =  s ./ (n-2)
    ss = ss / ((n-1)*(n-2))

    for (i, j, di, v) in DistanceMatrixDataIter(A)
        A.data[di] = i == j ? zero(T) : A.data[di] - s[i] - s[j] + ss
    end
    A.centering = :ucenter
    A.dvar = nothing
    A
end


"""
    dcov(A::DistanceMatrix{T}, B::DistanceMatrix{T}) where {T <: Real}

Computes the distance covariance of two centered DistanceMatrices `A` and `B`.
"""
function dcov(A::DistanceMatrix{T}, B::DistanceMatrix{T}) where {T <: Real}
    n = size(A, 1)
    n == size(B, 1) || throw(DimensionMismatch("matrices must be of same size"))
    A.centering != :none || throw("matrices must be centered")
    A.centering == B.centering || throw("matrices must be centered the same")

    s0 = zero(T)
    @simd for i in 1:n
        @inbounds s0 += A.data[i] * B.data[i]
    end
    s1 = zero(T)
    @simd for i in (n+1):length(A.data)
        @inbounds s1 += A.data[i] * B.data[i]
    end
    sqrt((s0 + 2 * s1) / (n*n))
end


"""
    dvar(A::DistanceMatrix{T}) where {T <: Real}

Computes the distance variance of a centered DistanceMatrices `A`.
Stores the variance alongside the DistanceMatrix for future use.
"""
function dvar(A::DistanceMatrix{T}) where {T <: Real}
    isnothing(A.dvar) || return(A.dvar)

    n = size(A, 1)
    A.centering != :none || throw("matrix must be centered")

    s0 = zero(T)
    @simd for i in 1:n
        @inbounds s0 += A.data[i] * A.data[i]
    end
    s1 = zero(T)
    @simd for i in (n+1):length(A.data)
        @inbounds s1 += A.data[i] * A.data[i]
    end
    A.dvar = sqrt((s0 + 2 * s1) / (n*n))
    A.dvar
end


"""
    dcor(A::DistanceMatrix{T}, B::DistanceMatrix{T})

Computes the distance correlation of two centered DistanceMatrices `A` and `B`.
"""
function dcor(A::DistanceMatrix{T}, B::DistanceMatrix{T}) where {T <: Real}
    dvarAdvarB = dvar(A) * dvar(B)
    iszero(dvarAdvarB) && return zero(T)
    dcov(A, B) / sqrt(dvarAdvarB)
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

0.3742040504583155
```
"""
function dcor(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Real}
    length(x) == length(y) || throw(DimensionMismatch("vectors must be of same length"))
    A = dcenter!(DistanceMatrix(x))
    B = dcenter!(DistanceMatrix(y))
    dcor(A, B)
end

end # module
