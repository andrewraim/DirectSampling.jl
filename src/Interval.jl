"""
This helper struct is used in the small_rects implementation of Stepdown.
It represents an interval [x,y] with function values h(x) and h(y).
"""
struct Interval
	log_x::Float64
	log_y::Float64
	log_hx::Float64
	log_hy::Float64
	function Interval(log_x, log_y, log_hx, log_hy)
		if log_x > log_y || log_hx < log_hy
			error("Require log_x <= log_y and log_hx >= log_hy for Interval")
		end
		new(convert(Float64, log_x),
			convert(Float64, log_y),
			convert(Float64, log_hx),
			convert(Float64, log_hy))
	end
end

function log_width(a::Interval)
	log_sub(a.log_y, a.log_x)
end

function log_height(a::Interval)
	log_sub(a.log_hx, a.log_hy)
end

function print(x::Interval; log_scale = false)
	if log_scale
		@printf "log_x: %g" x.log_x
		@printf "log_y: %g" x.log_y
		@printf "log_h_x: %g" x.log_hx
		@printf "log_h_y: %g" x.log_hy
		@printf "width: %g" log_width(x)
		@printf "height: %g" log_height(x)
	else
		@printf "x: %g\n" exp(x.log_x)
		@printf "y: %g\n" exp(x.log_y)
		@printf "h_x: %g\n" exp(x.log_hx)
		@printf "h_y: %g\n" exp(x.log_hy)
		@printf "width: %g\n" exp(log_width(x))
		@printf "height: %g\n" exp(log_height(x))
	end
end


