"""
Define abstract BaseDistribution

TBD: subclasses should implement the following functions
(DGeomBaseDistribution is an example of this type):

function logpdf(g::T, x) where {T<:DGeomBaseDistribution}
function pr_interval(g::T, x1, x2) where {T<:DGeomBaseDistribution}
function q_truncated(g::T, p, x_min, x_max) where {T<:DGeomBaseDistribution}
"""
abstract type BaseDistribution end

function pdf(g::T, x) where {T<:BaseDistribution}
	exp(log_pdf(g, x))
end

# Take one draw from this distribution after truncating to interval
# [x_min, x_max].
function r_truncated(g::T, x_min, x_max) where {T<:BaseDistribution}
	u = rand(Float64, 1)
	q_truncated(g, u[1], x_min, x_max)
end

function log_pdf(g::T, x) where {T<:BaseDistribution}
	error("log_pdf: please implement for your BaseDistribution")
end

function log_pr_interval(g::T, x1, x2) where {T<:BaseDistribution}
	error("log_pr_interval: please implement for your BaseDistribution")
end

function q_truncated(g::T, p, x_min, x_max) where {T<:BaseDistribution}
	error("q_truncated: please implement for your BaseDistribution")
end


