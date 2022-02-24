import DirectSampling: BaseDistribution, log_pdf, log_pr_interval, q_truncated
include("dgeom.jl")

"""
Define DGeomBaseDistribution
"""
struct DGeomBaseDistribution<:BaseDistribution
	rho::Float64
end

function DGeomBaseDistribution(rho)
	if rho < zero(rho) || rho > one(rho)
		error("rho must be between 0 and 1")
	end
	DGeomBaseDistribution(rho)
end

# Evaluate the DGeom(rho) density function
function log_pdf(g::T, x) where {T<:DGeomBaseDistribution}
	logpdf_dgeom(x, g.rho)
end

# Compute Pr(x1 < X <= x2) probability where X ~ DGeom(rho)
function log_pr_interval(g::T, x1, x2) where {T<:DGeomBaseDistribution}
	out = cdf_dgeom(x2, g.rho) - cdf_dgeom(x1, g.rho)
	log(out)
end

# Quantile function of DGeom truncated to [x_min, x_max]
function q_truncated(g::T, p, x_min, x_max) where {T<:DGeomBaseDistribution}
	p_min = cdf_dgeom(ceil(x_min) - 1, g.rho)
	p_max = cdf_dgeom(floor(x_max), g.rho)
	x = quantile_dgeom((p_max - p_min)*p + p_min, g.rho)
	max(ceil(x_min), min(x, floor(x_max)))
end

