"""
This class represents a step function approximation to the non-increasing
Pr(A_u) for u in [0,1].
"""
struct Stepdown
	w::WeightFunction
	g::BaseDistribution
	log_x_vals::Array{Float64}
	log_h_vals::Array{Float64}
	cum_probs::Array{Float64}
	norm_const::Float64
	max_jump::Float64
end

function log_p(w, g, log_u)
	endpoints = roots(w, log_c(w) + log_u)
	log_pr_interval(g, endpoints.first, endpoints.second)
end

function Stepdown(w, g; tol, N, method)
	# TBD: Any error checking??
	# if z < zero(z) || sigma2 < zero(sigma2)
	# 	error("z and sigma2 must be nonnegative")
	# end

	# TBD: let these be passed in!
	dist(x,y) = y - x
	midpoint(x,y) = (x + y) / 2

	# First, make sure the widest possible A_u intersects with [x_lo, x_hi]. If it
	# doesn't, the rest of the algorithm won't work, so bail out.
	log_prob_max = log_p(w, g, -Inf)
	if isinf(log_prob_max)
		error("Could not find any u such that P(X in A_u) > 0")
	end

	# Find L, the largest point where where P(A_L) = prob_max.
	# First find a finite lower bound for log_L.
	log_L_lo = -1
	log_L_hi = 0
	log_prob = log_p(w, g, log_L_lo)
	while log_prob < log_prob_max
		log_L_lo *= 2
		log_prob = log_p(w, g, log_L_lo)
	end

	# Do a bisection search to find log_L between log_L_lo and log_L_hi.

	pred1(log_u) = log_p(w, g, log_u) < log_prob_max
	delta = tol * (log_L_hi - log_L_lo)
	log_L = bisection(log_L_lo, log_L_hi, pred1, midpoint, dist, tol = delta)

	# Do a bisection search to find U, the smallest point where P(A_U) = 0.
	pred2(log_u) = isinf(log_p(w, g, log_u))
	delta = tol * (0 - log_L)
	log_U = bisection(log_L, 0, pred2, midpoint, dist, tol = delta)

	# Now fill in points between L and U
	if method == "equal_steps"
		out = equal_steps(w, g, exp(log_L), exp(log_U), log_prob_max, N)
		log_x_vals = out.first
		log_h_vals = out.second
	elseif method == "small_rects"
		out = small_rects(w, g, log_L, log_U, log_prob_max, N)
		log_x_vals = out.first
		log_h_vals = out.second
	else
		msg = string(
			@sprintf("Unknown method: %s. ", method),
			"Currently support equal_steps and small_rects\n")
		error(msg)
	end

	(areas, cum_probs, norm_const, max_jump) = initialize(log_x_vals, log_h_vals)
	Stepdown(w, g, log_x_vals, log_h_vals, cum_probs, norm_const, max_jump)
end

function equal_steps(w, g, L, U, log_prob_max, N)
	log_x_vals = Array{Float64}(undef, N+2)
	log_h_vals = Array{Float64}(undef, N+2)
	
	# Add pieces for the interval [-Inf, 0) where h(u) = 0, and
	# the interval [0, L) where h(u) = 1
	log_x_vals[1] = -Inf
	log_h_vals[1] = log_prob_max

	for i = 0:(N+1)
		prop = i / N
		log_u = log(L + prop * (U - L))
		log_x_vals[i+1] = log_u
		log_h_vals[i+1] = log_p(w, g, log_u)
	end

	return Pair(log_x_vals, log_h_vals)
end

function small_rects(w, g, log_L, log_U, log_prob_max, N)
	midpoint(x,y) = (x + y) / 2

	# This queue should be in max-heap order by height
	q = PriorityQueue{Interval,Float64}(Base.Order.Reverse)

	intvl = Interval(exp(log_L), exp(log_U), exp(log_p(w, g, log_L)), exp(log_p(w, g, log_U)))
	enqueue!(q, intvl, area(intvl))

	# Try to be efficient by preallocating log_x_vals and
	# log_h_vals to the (maximum) size needed, then copying these temporary
	# structures to _log_x_vals and _log_h_vals afterward. (It looks like
	# Rcpp NumericVectors cannot have their allocation controlled in this way).
	log_x_vals = Array{Float64}(undef, N+2)
	log_h_vals = Array{Float64}(undef, N+2)

	log_x_vals[1] = -Inf
	log_x_vals[2] = log_L
	log_x_vals[3] = log_U

	log_h_vals[1] = log_prob_max
	log_h_vals[2] = log_p(w, g, log_L)
	log_h_vals[3] = log_p(w, g, log_U)

	# We already have three (x, h(x)) pairs from above
	iter = 3

	while length(q) > 0 && iter < N+2
		# Get the interval with the largest priority
		int_top = dequeue!(q)

		# Break the interval int_top into two pieces: left and right.
		x_new = midpoint(int_top.x, int_top.y)
		log_x_new = log(x_new)
		log_h_new = log_p(w, g, log_x_new)
		h_new = exp(log_h_new)

		# Add the midpoint to our list of knots
		iter += 1
		log_x_vals[iter] = log_x_new
		log_h_vals[iter] = log_h_new

		# Add the interval which represents [int_top.x, x_new]
		int_left = Interval(int_top.x, x_new, int_top.fx, h_new)
		enqueue!(q, int_left, area(int_left))

		# Add the interval which represents [x_new, int_top.x]
		int_right = Interval(x_new, int_top.y, h_new, int_top.fy)
		enqueue!(q, int_right, area(int_right))
	end

	# Keep only the entries that we assigned and sort.
	# The x values should be ascending and the h values should be descending
	log_x_vals = sort(log_x_vals[1:iter])
	log_h_vals = sort(log_h_vals[1:iter], rev = true)

	return Pair(log_x_vals, log_h_vals)
end

function quantile(s::Stepdown, p)
	# Recall that cum_probs is a sorted vector
	if p < first(s.cum_probs)
		# p occurs before the first jump. There is no probability yet, so the
		# quantile is 0.
		out = 0
	elseif p >= last(s.cum_probs)
		out = 1
	else
		# p occurs after the first jump. Use find_interval to locate the two
		# cutpoints, then do a linear interpolation between them.
		j1 = find_interval(p, s.cum_probs)
		j2 = j1 + 1
		x1 = exp(s.log_x_vals[j1])
		x2 = exp(s.log_x_vals[j2])
		cp1 = s.cum_probs[j1]
		cp2 = s.cum_probs[j2]
		out = x1 + (x2 - x1) * (p - cp1) / (cp2 - cp1)

		if true
			@info "Searching for p = " p
			@info "s.cum_probs", s.cum_probs
			@info "j1 = ", j1
			@info "j2 = ", j2
			@info "x1 = ", x1
			@info "x2 = ", x2
			@info "exp(s.log_x_vals)", map(exp, s.log_x_vals)
			@info "cp1 = ", cp1
			@info "cp2 = ", cp2
			@info "out = ", out
		end
	end

	return out
end

function rand(s::Stepdown, n::Int)
	u = rand(Float64, n)
	map(x -> quantile(s,x), u)
end

function log_pdf(s::Stepdown, x; normalize = true)
	# Get the idx such that h_vals[idx] <= log(x) < h_vals[idx+1]
	idx = find_interval(log(x), s.log_x_vals)
	out = s.log_h_vals[idx]
	normalize ? out - s.norm_const : out
end

function pdf(s::Stepdown, x; normalize = true)
	exp(logpdf(s, x; normalize))
end


function cdf(s::Stepdown, x)
	if (log(x) > last(s.log_x_vals))
		out = 1
	else
		# Get the idx such that h_vals[idx] <= log(x) < h_vals[idx+1]	
		j1 = find_interval(log(x), s.log_x_vals)
		j2 = j1 + 1
		x1 = exp(s.log_x_vals[j1])
		x2 = exp(s.log_x_vals[j2])
		cp1 = s.cum_probs[j1]
		cp2 = s.cum_probs[j2]
		cp1 + (x - x1) / (x2 - x1) * (cp2 - cp1)
	end
end

# Update step function with new log_u and log_p values. Let's revisit
# implementing this one later. This is currently very inefficient, but
# we can consider making Stepdown mutable, and also using heaps to
# store log_x_vals and log_h_vals in sorted order.
function add(s::Stepdown, log_u)
	log_x_vals = s.log_x_vals
	log_h_vals = s.log_h_vals

	# Add log_u to s.log_x_vals. Ensure log_x_vals remains sorted in
	# nondecreasing order
	log_x_vals = sort(cat(s.log_x_vals, log_u, dims = 1))

	# Apply log_p to log_u and save the result to s.log_h_vals. Ensure
	# log_h_vals remains sorted in nonincreasing order	
	log_h = log_p(s.w, s.g, log_u)
	log_h_vals = sort(cat(s.log_h_vals, log_h, dims = 1), rev = true)

	(areas, cum_probs, norm_const, max_jump) = initialize(log_x_vals, log_h_vals)
	Stepdown(s.w, s.g, log_x_vals, log_h_vals, cum_probs, norm_const, max_jump)
end

# Compute normalizing constant, cumulative probabilities, and maximum jump
# based on log_x_vals and log_h_vals.
function initialize(log_x_vals, log_h_vals)
	N = length(log_x_vals) - 2
	areas = Vector{Float64}(undef, N+1)

	cum_probs = Vector{Float64}(undef, N+1)
	norm_const = 0
	max_jump = 0

	if false
		# This is a more readable version of the calculation that may be less precise
		for i = 1:(N+1)
			h_cur = exp(log_h_vals[i])
			h_next = exp(log_h_vals[i+1])
			x_cur = exp(log_x_vals[i])
			x_next = exp(log_x_vals[i+1])
			areas[i] = h_cur * (x_next - x_cur)
			norm_const += areas[i]
			cum_probs[i] = norm_const

			# Avoid division by zero, which can occur in later exp(h) values
			max_jump = max(max_jump, (h_next > 0) * h_cur/h_next)
		end
	else
		# This version tries to be more precise and do more work on the log-scale
		for i = 1:(N+1)
			log_h_cur = log_h_vals[i]
			log_h_next = log_h_vals[i+1]
			log_x_cur = log_x_vals[i]
			log_x_next = log_x_vals[i+1]
			areas[i] = exp(log_x_next + log_h_cur) - exp(log_x_cur + log_h_cur)
			norm_const += areas[i]
			cum_probs[i] = norm_const

			# Avoid division by zero, which can occur in later exp(h) values
			max_jump = max(max_jump, !isinf(log_h_next) * exp(log_h_cur - log_h_next))
		end
	end

	cum_probs = cum_probs / norm_const;

	if any(map(isnan, cum_probs))
		error("NA values found in cum_probs. Approximation may have failed")
	end

	return areas, cum_probs, norm_const, max_jump
end

