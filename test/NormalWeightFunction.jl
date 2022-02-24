using DirectSampling
import DirectSampling.WeightFunction

"""
Define NormalWeightFunction
"""
struct NormalWeightFunction <: WeightFunction
	z::Float64
	mu::Float64
	sigma2::Float64
end

# Return the roots of the equation w(x) = a, which is equivalent to
# log(w(x)) = log(a). Roots are returned in increasing order.
function roots(w::T, log_a) where {T<:NormalWeightFunction}
	if -2 * w.sigma2 * log_a < 0 && log_a <= log_c(w)
		x1 = w.z - w.mu
		x2 = w.z - w.mu
	else
		x1 = w.z - (w.mu + sqrt(-2 * w.sigma2 * log_a))
		x2 = w.z - (w.mu - sqrt(-2 * w.sigma2 * log_a))
	end
	Pair(x1, x2)
end

function log_c(w::T) where {T<:NormalWeightFunction}
	return 0
end

function log_eval(w::T, x) where {T<:NormalWeightFunction}
	-1 / (2*w.sigma2) * (w.z - x - w.mu)^2
end


