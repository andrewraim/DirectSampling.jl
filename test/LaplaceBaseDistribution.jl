import DirectSampling: BaseDistribution, log_pdf, log_pr_interval, q_truncated
include("laplace.jl")

"""
Define LaplaceBaseDistribution
"""
struct LaplaceBaseDistribution<:BaseDistribution
	lambda::Float64
end

function LaplaceBaseDistribution(lambda)
	if lambda < zero(lambda)
		error("lambda must be nonnegative")
	end
	LaplaceBaseDistribution(lambda)
end

# Evaluate the Laplace(0, lambda) density function
function log_pdf(g::T, x) where {T<:LaplaceBaseDistribution}
	logpdf_laplace(x, 0, g.lambda)
end

# Compute Pr(x1 < X <= x2) probability where X ~ Laplace(0, lambda)
function log_pr_interval(g::T, x1, x2) where {T<:LaplaceBaseDistribution}
	out = cdf_laplace(x2, 0, g.lambda) - cdf_laplace(x1, 0, g.lambda)
	log(out)
end

# Quantile function of Laplace(0, lambda) truncated to [x_min, x_max]
function q_truncated(g::T, p, x_min, x_max) where {T<:LaplaceBaseDistribution}
	p_min = cdf_laplace(x_min, 0, g.lambda)
	p_max = cdf_laplace(x_max, 0, g.lambda)
	x = quantile_laplace((p_max - p_min)*p + p_min, 0, g.lambda);
	max(x_min, min(x, x_max));
end

