"""
    interp1d(x::Number, x1::Number, x2::Number, y1, y2)

Approximate `y(x)` with a linear interpolation between two points `(x1, y1)` and `(x2, y2)`. This method broadcasts over `y1` and `y2` so they can be vector-valued.

# Examples
```jldoctest; setup = :(import BolometricCorrections: interp1d)
julia> interp1d(1.5, 1.0, 2.0, 1:10, 2:11) ≈ 1.5:10.5
true
```
"""
interp1d(x::Number, x1::Number, x2::Number, y1, y2) = @. (y1 * (x2 - x) + y2 * (x - x1)) / (x2 - x1)
"""
    interp2d(x::Number, y::Number, x1::Number,
             x2::Number, y1::Number, y2::Number, z1_1, z2_1, z1_2, z2_2)

Approximate `z(x,y)` with a linear interpolation between the four points

 - `(x1, y1, z1_1)`
 - `(x2, y1, z2_1)`
 - `(x1, y2, z1_2)`
 - `(x2, y2, z2_2)`

This method broadcasts over `y1` and `y2` so they can be vector-valued.

# Examples
```jldoctest; setup = :(import BolometricCorrections: interp2d)
julia> interp2d(1.75, 1.5, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 3.0, 4.0) ≈ 2.75
true
```
"""
function interp2d(x::Number, y::Number, x1::Number, x2::Number, y1::Number, y2::Number, z1_1, z2_1, z1_2, z2_2)
    return @. (z1_1 * (x2 - x) * (y2 - y) + z1_2 * (x2 - x) * (y - y1) + z2_1 * (x - x1) * (y2 - y) + z2_2 * (x - x1) * (y - y1)) / (x2 - x1) / (y2 - y1)
end

# Use statically known size from filters argument to repack submatrix
# into a vector of SVectors to pass into interpolator
"""
repack_submatrix(submatrix::AbstractArray{T},
                 dim1::Int, dim2::Int,
                 filters::Union{NTuple{N}, Val{N}}) where {T, N}
Use static size information (`N`) to reshape `submatrix` for interpolation.
First, reshape `submatrix` which has size `(dim1 * dim2, N)` to `(dim1, dim2, N)`. Then 
convert to a `Vector{SVector{N, T}}` for more efficient use with Interpolations.jl.
"""
function repack_submatrix(submatrix::AbstractArray, dim1::Int, dim2::Int, ::NTuple{N}) where {N}
    return repack_submatrix(submatrix, dim1, dim2, Val(N))
end
function repack_submatrix(submatrix::AbstractArray{T}, dim1::Int, dim2::Int, ::Val{N}) where {T, N}
    submatrix = reshape(submatrix,
                        dim1, # length(logg),
                        dim2, # length(teff),
                        N)
    return [SVector{N, T}(view(submatrix,i,j,:)) for i=axes(submatrix,1), j=axes(submatrix,2)]
end