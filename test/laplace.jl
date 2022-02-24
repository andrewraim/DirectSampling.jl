using Distributions

"""
It might be better to use object-orientation like in Distributions.jl
but let's start with something simpler. We don't need vectorized
versions of each function, except in cases where calling 'map' would
be very inefficient.
"""

function insupport_laplace(x)
	typeof(x) <: Real
end

function logpdf_laplace(x, mu, lambda)
	-log(2) - log(lambda) - abs(x - mu) / lambda
end

function pdf_laplace(x, mu, lambda)
	exp(logpdf_laplace(x, mu, lambda))
end

function cdf_laplace(x, mu, lambda)
	if x <= mu
		0.5 * exp((x - mu) / lambda)
	else
		1 - 0.5 * exp(-(x - mu) / lambda)
	end
end

function rand_laplace(mu, lambda)
	d = Laplace(mu, lambda)
	x = rand(d,1) - rand(d,1)
	x[1]
end

function rand_laplace(n, mu, lambda)
	d = Laplace(mu, lambda)
	rand(d,n) - rand(d,n)
end

function quantile_laplace(q, mu, lambda)
	if q <= 0.5
		mu + lambda * log(2*q)
	else
		mu - lambda * log(2 - 2*q)
	end
end


