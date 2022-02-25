function direct_sampler_ar(n, w::T1, g::T2; tol, N, max_rejections, fill_method,
priority_weight) where {T1<:WeightFunction, T2<:BaseDistribution}
	# Use our Stepdown approximation to draw from p(u)
	step = Stepdown(w, g; tol = tol, N = N, method = fill_method,
		priority_weight = priority_weight)

	log_u = Array{Real}(undef, n)
	rejections = zero(Int64)
	accept = false

	# Because the step function is an upper bound for P(A_u), the constant M in
	# the acceptance ratio is always M = 1.
	log_M = 0

	for i = 1:n
		accept = false
		while !accept && rejections < max_rejections
			v = rand(Float64, 1)
			log_u_proposal = rand(step, n = 1, take_log = true)
			log_p_val = log_p(step, log_u_proposal)
			log_h_val = pdf(step, log_u_proposal, true, take_log = false)
			log_ratio = log_p_val - log_h_val - log_M

			if log(v) < log_ratio
				# Accept u as a draw from p(u)
				log_u[i] = log_u_proposal
				accept = true
			else
				# Reject u and add it to knots
				add(step, log_u_proposal)
				rejections += 1
			end
		end
	end

	if rejections == max_rejections
		msg = @sprintf "Reached maximum number of rejections: %d" max_rejections
		error(msg)
	end

	# Draw from g(x | u) for each value of log(u)
	for i = 1:n
		endpoints = roots(w, log_c(w) + log_u[i])
		x[i] = r_truncated(g, endpoints.first, endpoints.second)
	end

	return x, rejections
end

