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
# @. (mat1 * (feh2 - feh) + mat2 * (feh - feh1)) / (feh2 - feh1)
