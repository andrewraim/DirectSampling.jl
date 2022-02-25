function direct_sampler(w::T1, g::T2; n = 1, tol = 1e-10, N = 50,
knot_method = SMALL_RECTS, priority_weight = 1/2) where {T1<:WeightFunction,
T2<:BaseDistribution}

	# For now, let's assume geometric midpoint
	log_midpoint(log_x, log_y) = geometric_midpoint(log_x, log_y; take_log = true)

	# Use our Stepdown approximation to draw from p(u)
	step = Stepdown(w, g; tol = tol, N = N, knot_method = knot_method,
		priority_weight = priority_weight, log_midpoint = log_midpoint)
	log_u = rand(step; n = n, take_log = true)

	# Draw from g(x | u) for each value of log(u)
	x = zeros(Float64, n)
	for i = 1:n
		endpoints = roots(w, log_c(w) + log_u[i])
		x[i] = r_truncated(g, endpoints.first, endpoints.second)
	end

	return x
end

