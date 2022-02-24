"""
Assume cutpoints is a sorted vector, where elements represent endpoints
of adjacent intervals [c_1,c_2), [c_2,c_3) ..., [c_{N-1},c_{N}). Return:
- the index i such that z is in [c_i,c_{i+1}), or
- return 0 if z < c_{1}, or
- N if z > c_{N}.

TBD: Should we take away any of the type specifications?

julia> cc = collect(1:3)
julia> print(cc)
julia> find_interval(2.5, cc)
"""
function find_interval(z::S, cutpoints::Array{T,1}) where {S<:Real, T<:Real}
	N = length(cutpoints)

	if z < first(cutpoints)
		return 0
	elseif z >= last(cutpoints)
		return N
	end

	find_interval_pred(x) = cutpoints[Int(floor(x))] >= z
	floor_midpoint(x,y) = floor((y + x) / 2)
	interval_length(x,y) = y - x

	out = bisection(1, N, find_interval_pred, floor_midpoint, interval_length; tol = 1)
	return Int(floor(out))
end

