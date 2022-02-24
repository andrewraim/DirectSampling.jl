"""
This helper struct is used in the small_rects implementation of Stepdown.
It represents an interval [x,y] with function values h(x) and h(y).
"""
struct Interval
	x::Float64
	y::Float64
	fx::Float64
	fy::Float64
	function Interval(x, y, fx, fy)
		if x > y || fx < fy
			error("x <= y and fx >= fy must be true for a Stepdown Interval")
		end
		new(convert(Float64, x),
			convert(Float64, y),
			convert(Float64, fx),
			convert(Float64, fy))
	end
end

function area(a::Interval)
	(a.y - a.x) * (a.fx - a.fy)
end

