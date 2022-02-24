using Distributions

"""
It might be better to use object-orientation like in Distributions.jl
but let's start with something simpler. We don't need vectorized
versions of each function, except in cases where calling 'map' would
be very inefficient.
"""

function insupport_dgeom(x)
	isinteger(x)
end

function logpdf_dgeom(x, p)
	if insupport_dgeom(x)
		log(p) - log(2-p) + abs(x) * log(1-p)
	else
		log(zero(p))
	end
end

function pdf_dgeom(x, p)
	exp(logpdf_dgeom(x, p))
end

function cdf_dgeom(x, p)
	xx = floor(x)

	if xx >= 0
		(1-p) / (2-p) + (1 - (1-p)^(xx+1)) / (2-p)
	else
		(1-p)^(-xx) / (2-p)
	end
end

function rand_dgeom(p)
	d = Geometric(p)
	x = rand(d,1) - rand(d,1)
	x[1]
end

function rand_dgeom(n, p)
	d = Geometric(p)
	rand(d,n) - rand(d,n)
end

function quantile_dgeom(q, p)
	x_neg = ceil( -log((2-p)*q) / log(1-p) )
	x_pos = ceil( (log(1-q) + log(2-p)) / log(1-p) - 1 )
	x_neg*(x_neg < 0) + x_pos*(x_pos >= 0)
end


