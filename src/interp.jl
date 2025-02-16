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
