function direct_sampler(n, w::T1, g::T2; tol, N, fill_method) where
{T1<:WeightFunction, T2<:BaseDistribution}
	# Use our Stepdown approximation to draw from p(u)
	step = Stepdown(w, g; tol = tol, N = N, method = fill_method,
		priority_weight = priority_weight)
	log_u = rand(step; n = n, take_log = true)

	# Draw from g(x | u) for each value of log(u)
	x = Array{Real}(undef, n)
	for i = 1:n
		endpoints = roots(w, log_c(w) + log_u[i])
		x[i] = r_truncated(g, endpoints.first, endpoints.second)
	end

	return x
end

