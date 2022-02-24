# using DirectSampling
import DirectSampling: WeightFunction, log_c, log_eval, roots

"""
Define LognormalWeightFunction
"""
struct LognormalWeightFunction <: WeightFunction
	z::Float64
	mu::Float64
	sigma2::Float64
end

function LognormalWeightFunction(z, mu, sigma2)
	if sigma2 < zero(sigma2)
		error("sigma2 must be nonnegative")
	end
	LognormalWeightFunction(
		convert(Float64, z),
		convert(Float64, mu),
		convert(Float64, sigma2))
end


# Return the roots of the equation w(x) = a, which is equivalent to
# log(w(x)) = log(a). Roots are returned in increasing order.
function roots(w::T, log_a) where {T<:LognormalWeightFunction}
	if w.sigma2 < 2*(w.mu + log_a) && log_a <= log_c(w)
		x1 = w.z - exp(w.mu - w.sigma2)
		x2 = w.z - exp(w.mu - w.sigma2)
	else
		x1 = w.z - exp((w.mu - w.sigma2) + sqrt(w.sigma2 * (w.sigma2 - 2*(w.mu + log_a))))
		x2 = w.z - exp((w.mu - w.sigma2) - sqrt(w.sigma2 * (w.sigma2 - 2*(w.mu + log_a))))
	end

	# Edge case: If z is an integer and the larger root is numerically close to
	# z, make the root smaller by an epsilon. This helps to ensure that the
	# endpoint z is not part of the support of p(u).
	if w.z - floor(w.z) < 1e-10 && w.z - x2 < 1e-6
		x2 = w.z - 1e-6
	end
	x1 = min(x1, x2)

	Pair(x1, x2)
end

function log_c(w::T) where {T<:LognormalWeightFunction}
	-(w.mu - w.sigma2) / 2
end

function log_eval(w::T, x) where {T<:LognormalWeightFunction}
	out = -Inf
	if x < w.z
		out = -log(w.z - x) - (log(w.z - x) - w.mu)^2 / (2 * w.sigma2)
	end
	return out
end

