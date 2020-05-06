var documenterSearchIndex = {"docs":
[{"location":"#EnergyStatistics.jl-1","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"","category":"section"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"In statistics distance correlation or distance covariance is a measure of dependence between two paired random vectors. The population distance correlation coefficient is zero if and only if the random vectors are independent. Thus, distance correlation measures both linear and nonlinear association between two random vectors. This is in contrast to Pearson's correlation, which can only detect linear association between two random variables. See here for references and more details.","category":"page"},{"location":"#Installation-1","page":"EnergyStatistics.jl","title":"Installation","text":"","category":"section"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"This package can be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run","category":"page"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"pkg> add EnergyStatistics","category":"page"},{"location":"#General-Usage-1","page":"EnergyStatistics.jl","title":"General Usage","text":"","category":"section"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"Given two vectors x and y the distance correlation dcor can simply computed:","category":"page"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"using EnergyStatistics\n\nx = collect(-1:0.01:1)\ny = map(x -> x^4 - x^2, x)\n\ndcor(x, y) ≈ 0.374204050","category":"page"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"These two vectors are clearly associated. However, their (Pearson) correlation coefficient vanishes suggesting that they are independent. The finite distance correlation dcor reveals their non-linear association.","category":"page"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"Function to compute the distance covariance dcov and distance variance dvar are also supplied.","category":"page"},{"location":"#Advanced-Usage-1","page":"EnergyStatistics.jl","title":"Advanced Usage","text":"","category":"section"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"The computation of a 'DistanceMatrix' is computationally expansive. Especially the computation of n(n-1)/2 pairwise distances for vectors of length n and the subsequent centering of the distance matrix take time and memory. In cases where one wants to compute several distance correlations and keep intermediate results of the distance computations and centering one can do so. For example:","category":"page"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"Dx = EnergyStatistics.dcenter!(EnergyStatistics.DistanceMatrix(x))\nDy = EnergyStatistics.dcenter!(EnergyStatistics.DistanceMatrix(y))\nDz = EnergyStatistics.dcenter!(EnergyStatistics.DistanceMatrix(z))\n\ndcor_xy = dcor(Dx, Dy)\ndcor_xz = dcor(Dx, Dz)","category":"page"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"will run faster than","category":"page"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"dcor_xy = dcor(x, y)\ndcor_xz = dcor(x, z)","category":"page"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"since the distance matrix Dx for the vector x is only computed once.","category":"page"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"You can also construct distance matrices using other distance measures than the (default) abs.","category":"page"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"AA = EnergyStatistics.dcenter!(EnergyStatistics.DistanceMatrix(Float64, x, abs2))","category":"page"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"Instead of double centering via dcenter! one may also use U-centering via the ucenter! function.","category":"page"},{"location":"#References-1","page":"EnergyStatistics.jl","title":"References","text":"","category":"section"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"See the wikipedia page for references and more details.","category":"page"},{"location":"#Functions-1","page":"EnergyStatistics.jl","title":"Functions","text":"","category":"section"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"EnergyStatistics.dcor(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Real}","category":"page"},{"location":"#EnergyStatistics.dcor-Union{Tuple{T}, Tuple{AbstractArray{T,1},AbstractArray{T,1}}} where T<:Real","page":"EnergyStatistics.jl","title":"EnergyStatistics.dcor","text":"dcor(x::AbstractVector{T}, y::AbstractVector{T}) where T <: Real\n\nComputes the distance correlation of samples x and y.\n\nusing EnergyStatistics\nx = collect(-1:0.01:1)\ny = @. x^4 - x^2\ndcor(x, y)\n\n# output\n\n0.3742040504583155\n\n\n\n\n\n","category":"method"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"EnergyStatistics.dcov(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Real}","category":"page"},{"location":"#EnergyStatistics.dcov-Union{Tuple{T}, Tuple{AbstractArray{T,1},AbstractArray{T,1}}} where T<:Real","page":"EnergyStatistics.jl","title":"EnergyStatistics.dcov","text":"dcov(x::AbstractVector{T}, y::AbstractVector{T}) where T <: Real\n\nComputes the distance covariance of samples x and y.\n\n\n\n\n\n","category":"method"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"EnergyStatistics.dvar(x::AbstractVector{T}) where {T <: Real}","category":"page"},{"location":"#EnergyStatistics.dvar-Union{Tuple{AbstractArray{T,1}}, Tuple{T}} where T<:Real","page":"EnergyStatistics.jl","title":"EnergyStatistics.dvar","text":"dvar(x::AbstractVector{T}) where T <: Real\n\nComputes the distance variance of a sample x.\n\n\n\n\n\n","category":"method"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"EnergyStatistics.DistanceMatrix(x::AbstractVector{T}, dist = abs ) where {T}","category":"page"},{"location":"#EnergyStatistics.DistanceMatrix-Union{Tuple{AbstractArray{T,1}}, Tuple{T}, Tuple{AbstractArray{T,1},Any}} where T","page":"EnergyStatistics.jl","title":"EnergyStatistics.DistanceMatrix","text":"DistanceMatrix(x::AbstractVector{T}, dist = abs) where {T}\n\nComputes the matrix of pairwise distance of x. The distance measure dist is abs as default.\n\nusing EnergyStatistics\nx = [1.0, 2.0]\nEnergyStatistics.DistanceMatrix(x)\n\n# output\n\n2×2 EnergyStatistics.DistanceMatrix{Float64}:\n 0.0  1.0\n 1.0  0.0\n\n\n\n\n\n","category":"method"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"EnergyStatistics.dcenter!(A::EnergyStatistics.DistanceMatrix{T}) where {T <: Real}","category":"page"},{"location":"#EnergyStatistics.dcenter!-Union{Tuple{EnergyStatistics.DistanceMatrix{T}}, Tuple{T}} where T<:Real","page":"EnergyStatistics.jl","title":"EnergyStatistics.dcenter!","text":"dcenter!(A::DistanceMatrix{T}) where {T <: Real}\n\nComputes the double centered matrix of A in place.\n\n\n\n\n\n","category":"method"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"EnergyStatistics.ucenter!(A::EnergyStatistics.DistanceMatrix{T}) where {T <: Real}","category":"page"},{"location":"#EnergyStatistics.ucenter!-Union{Tuple{EnergyStatistics.DistanceMatrix{T}}, Tuple{T}} where T<:Real","page":"EnergyStatistics.jl","title":"EnergyStatistics.ucenter!","text":"ucenter!(A::DistanceMatrix{T}) where {T <: Real}\n\nComputes the u-centered matrix of A in place.\n\n\n\n\n\n","category":"method"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"EnergyStatistics.dcor(A::EnergyStatistics.DistanceMatrix{T}, B::EnergyStatistics.DistanceMatrix{T}) where {T <: Real}","category":"page"},{"location":"#EnergyStatistics.dcor-Union{Tuple{T}, Tuple{EnergyStatistics.DistanceMatrix{T},EnergyStatistics.DistanceMatrix{T}}} where T<:Real","page":"EnergyStatistics.jl","title":"EnergyStatistics.dcor","text":"dcor(A::DistanceMatrix{T}, B::DistanceMatrix{T})\n\nComputes the distance correlation of two centered DistanceMatrices A and B.\n\n\n\n\n\n","category":"method"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"EnergyStatistics.dcov(A::EnergyStatistics.DistanceMatrix{T}, B::EnergyStatistics.DistanceMatrix{T}) where {T <: Real}","category":"page"},{"location":"#EnergyStatistics.dcov-Union{Tuple{T}, Tuple{EnergyStatistics.DistanceMatrix{T},EnergyStatistics.DistanceMatrix{T}}} where T<:Real","page":"EnergyStatistics.jl","title":"EnergyStatistics.dcov","text":"dcov(A::DistanceMatrix{T}, B::DistanceMatrix{T}) where {T <: Real}\n\nComputes the distance covariance of two centered DistanceMatrices A and B.\n\n\n\n\n\n","category":"method"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"EnergyStatistics.dvar(A::EnergyStatistics.DistanceMatrix{T}) where {T <: Real}","category":"page"},{"location":"#EnergyStatistics.dvar-Union{Tuple{EnergyStatistics.DistanceMatrix{T}}, Tuple{T}} where T<:Real","page":"EnergyStatistics.jl","title":"EnergyStatistics.dvar","text":"dvar(A::DistanceMatrix{T}) where {T <: Real}\n\nComputes the distance variance of a centered DistanceMatrices A. Stores the variance alongside the DistanceMatrix for future use.\n\n\n\n\n\n","category":"method"},{"location":"#Index-1","page":"EnergyStatistics.jl","title":"Index","text":"","category":"section"},{"location":"#","page":"EnergyStatistics.jl","title":"EnergyStatistics.jl","text":"","category":"page"}]
}