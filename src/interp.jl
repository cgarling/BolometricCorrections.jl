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

"""
    fill_bad_values(mat::AbstractMatrix{T}; 
                    isbad = Base.Fix1(==, zero(T)), 
                    window::Int = 1,
                    diag::Bool = false) where {T <: Number}
Returns a matrix similar to `mat` but with elements where `isbad(mat[i,j]) == true` are replaced by an average of valid neighbors. If `diag == true`, then elements on diagonals are included when calculating the replacement values; if `false`, then only elements in cardinal directions are considered. 
"""
function fill_bad_values(mat::AbstractMatrix{T}; 
                         isbad = Base.Fix1(==, zero(T)), 
                         window::Int = 1,
                         diag::Bool = false) where {T <: Number}
    a1, a2 = axes(mat)
    newmat = copy(mat)

    if diag # If counting elements in diagonal directions,
        steps = [(i, j) for i=-window:window, j=-window:window if (i ≠ 0) || (j ≠ 0)]
    else # If using only cardinal directions,
        steps = [(w, 0) for w in -window:window if w ≠ 0] ∪ [(0, w) for w in -window:window if w ≠ 0]
    end

    for i=a1, j=a2
        if isbad(mat[i, j])
            vals = T[]
            # for di in -window:window, dj in -window:window
            for (di, dj) in steps
                ii, jj = i + di, j + dj
                if checkbounds(Bool, mat, ii, jj) && !isbad(mat[ii, jj])
                    push!(vals, mat[ii, jj])
                end
            end
            if !isempty(vals)
                newmat[i, j] = mean(vals)
            end
        end
    end

    # Recursively call until no bad values remaining
    while any(isbad, newmat)
        newmat = fill_bad_values(newmat)
    end

    return newmat

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
    submatrix = reshape(submatrix, dim1, dim2, N)
    return [SVector{N, T}(view(submatrix,i,j,:)) for i=axes(submatrix,1), j=axes(submatrix,2)]
end

"""
repack_submatrix(submatrix::AbstractArray{T},
                 filters::Union{NTuple{N}, Val{N}}) where {T, N}
Use static size information (`N`) to reshape `submatrix` for interpolation.
The last dimension (which must have size `N`) is lifted out of the `submatrix`,
returning an array with one fewer dimension that contains elements of `SVector{N, T}}`
for more efficient use with Interpolations.jl.

`submatrix` must already be shaped to have the proper dimensionality for interpolation
after removal of the last dimension.
"""
repack_submatrix(submatrix::AbstractArray) = repack_submatrix(submatrix, Val(size(submatrix)[end]))
function repack_submatrix(submatrix::AbstractArray{T}, ::Val{N}) where {T, N}
    @argcheck size(submatrix)[end] == N "Last dimension must be $N"
    leading_dims = size(submatrix)[begin:end-1]
    out = Array{SVector{N, T}}(undef, leading_dims...)
    for I in CartesianIndices(out)
        out[I] = SVector{N, T}(view(submatrix, I, :)) # out[I] = SVector{N, T}(view(submatrix, Tuple(I)..., :))
    end
    return out
end