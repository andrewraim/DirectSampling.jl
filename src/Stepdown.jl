"""
This class represents a step function approximation to the non-increasing
Pr(A_u) for u in [0,1].
"""
@enum StepdownKnotMethod begin
	EQUAL_STEPS
	SMALL_RECTS
end

mutable struct Stepdown
	w::WeightFunction
	g::BaseDistribution
	log_x_vals::Array{Float64}
	log_h_vals::Array{Float64}
	cum_probs::Array{Float64}
	norm_const::Float64
	N::Integer
	tol::Real
	knot_order::Array{Integer}
	priority_weight::Real
	knot_method::StepdownKnotMethod
	log_midpoint
end

function log_p(w, g, log_u)
	endpoints = roots(w, log_c(w) + log_u)
	log_pr_interval(g, endpoints.first, endpoints.second)
end

function rects(s; take_log = false)
	k = s.N + 2
	log_x1 = s.log_x_vals[-k]
	log_x2 = s.log_x_vals[-1]
	log_h1 = s.log_h_vals[-k]
	log_h2 = s.log_h_vals[-1]
	log_area = zeros(k-1)

	for i in 1:(k-1)
		log_area[i] = log_sub(log_x2[i], log_x1[i]) + log_sub(log_h1[i], log_h2[i])
	end

	if take_log
		rects = Dict(
			"log_x1" => log_x1,
			"log_x2" => log_x2,
			"log_h1" => log_h1,
			"log_h2" => log_h2,
			"log_area" => log_area
		)
	else
		rects = Dict(
			"x1" => exp(log_x1),
			"x2" => exp(log_x2),
			"h1" => exp(log_h1),
			"h2" => exp(log_h2),
			"area" => exp(log_area)
		)
	end

	return rects
end

function Stepdown(w, g; tol = 1e-10, N = 50, knot_method = SMALL_RECTS,
priority_weight = 0.5,
log_midpoint = geometric_midpoint(log_x, log_y; take_log = true))
	# TBD: Any error checking??
	# if z < zero(z) || sigma2 < zero(sigma2)
	# 	error("z and sigma2 must be nonnegative")
	# end

	log_distance(log_x, log_y) = univariate_distance(log_x, log_y, take_log = true)

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

	pred_logL(log_u) = log_p(w, g, log_u) < log_prob_max
	log_delta1 = min(log_L_lo, log(tol))
	log_L = bisection(log_L_lo, log_L_hi, pred_logL, log_midpoint,
		log_distance, tol = log_delta1)

	# Do a bisection search to find U, the smallest point where P(A_U) = 0.
	pred_logU(log_u) = log_p(w, g, log_u) < log(1e-10) + log_p(w, g, log_L)
	log_delta2 = min(log_L, log(tol))
	log_U = bisection(log_L, 0, pred_logU, log_midpoint, log_distance,
		tol = log_delta2)

	# Now fill in points between L and U
	if knot_method == EQUAL_STEPS
		(a,b,c) = equal_steps(w, g, log_L, log_U, log_prob_max, N)
		log_x_vals = a
		log_h_vals = b
		knot_order = c
	elseif knot_method == SMALL_RECTS
		(a,b,c) = small_rects(w, g, log_L, log_U, log_prob_max,
			N = N, tol = tol, priority_weight = priority_weight,
			log_midpoint = log_midpoint)
		log_x_vals = a
		log_h_vals = b
		knot_order = c
	else
		error("knot_method must be equal_steps or small_rects")
	end

	(areas, cum_probs, norm_const) = steps2probs(log_x_vals, log_h_vals)
	return Stepdown(w, g, log_x_vals, log_h_vals, cum_probs, norm_const, N,
		tol, knot_order, priority_weight, knot_method, log_midpoint)
end

function equal_steps(w, g, log_L, log_U, log_prob_max, N)
	log_x_vals = Array{Float64}(undef, N+2)
	log_h_vals = Array{Float64}(undef, N+2)
	knot_order = collect(1:(N+2))

	# Add pieces for the interval [-Inf, 0) where h(u) = 0, and
	# the interval [0, L) where h(u) = 1
	log_x_vals[1] = -Inf
	log_h_vals[1] = log_prob_max

	# Compute log(U - L)
	log_span = log_sub(log_U, log_L)

	for i = 0:N
		log_prop = log(i / N)
		# Compute log_u = log(L + prop * (U - L))
		log_u = log_add(log_prop + log_span, log_L)
		log_x_vals[i+1] = log_u
		log_h_vals[i+1] = log_p(w, g, log_u)
	end

	return log_x_vals, log_h_vals, knot_order
end

function small_rects(w, g, log_L, log_U, log_prob_max; N, tol, priority_weight,
	log_midpoint)

	pw = priority_weight

	# This queue should be in max-heap order by height
	q = PriorityQueue{Interval,Float64}(Base.Order.Reverse)

	intvl = Interval(log_L, log_U, log_p(w, g, log_L), log_p(w, g, log_U))
	priority = pw * log_height(intvl) + (1-pw) * log_width(intvl)
	enqueue!(q, intvl, priority)

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
		log_x_new = log_midpoint(int_top.log_x, int_top.log_y)
		log_h_new = log_p(w, g, log_x_new)

		# Add the midpoint to our list of knots
		iter += 1
		log_x_vals[iter] = log_x_new
		log_h_vals[iter] = log_h_new

		# Add the interval which represents [int_top.x, x_new]
		int_left = Interval(int_top.log_x, log_x_new, int_top.log_hx, log_h_new)
		priority = pw * log_height(int_left)+ (1-pw) * log_width(int_left)
		enqueue!(q, int_left, priority)

		# Add the interval which represents [x_new, int_top.x]
		int_right = Interval(log_x_new, int_top.log_y, log_h_new, int_top.log_hy)
		priority = pw * log_height(int_right)+ (1-pw) * log_width(int_right)
		enqueue!(q, int_right, priority)
	end

	# Keep only the entries that we assigned and sort.
	# The x values should be ascending and the h values should be descending
	idx = sortperm(log_x_vals)
	log_x_vals_srt = log_x_vals[idx]
	log_h_vals_srt = log_h_vals[idx]
	knot_order = idx

	return log_x_vals_srt, log_h_vals_srt, knot_order
end

function quantile(s::Stepdown, p; take_log = false)

	# @printf "In quantile for Stepdown, p = %g and take_log = %d" p take_log
	# @printf "\n"

	cum_probs_ext = cat(0, s.cum_probs, dims = (1,1))
	idx = find_interval(p, cum_probs_ext)

	# Recall that cum_probs is a sorted vector
	if idx == 0
		# p occurs before the first jump. There is no probability yet, so the
		# quantile is 0.
		out = -Inf
	elseif idx >= s.N + 2
		out = 0
	else
		# p occurs after the first jump. Use find_interval to locate the two
		# cutpoints, then do a linear interpolation between them.
		j1 = idx
		j2 = idx + 1
		cp1 = cum_probs_ext[j1]
		cp2 = cum_probs_ext[j2]
		log_x1 = s.log_x_vals[j1]
		log_x2 = s.log_x_vals[j2]

		# Do the following computation, but be careful to keep x values on
		# the log-scale, since they may be extremely small.
		# out = x1 + (x2 - x1) * (p - cp1) / (cp2 - cp1)
		log_ratio = log(p - cp1) - log(cp2 - cp1)
		out = log_add(log_sub(log_x2, log_x1) + log_ratio, log_x1)

		# if true
		#	@info "Searching for p = " p
		#	@info "s.cum_probs", s.cum_probs
		#	@info "j1 = ", j1
		#	@info "j2 = ", j2
		#	@info "x1 = ", x1
		#	@info "x2 = ", x2
		#	@info "exp(s.log_x_vals)", map(exp, s.log_x_vals)
		#	@info "cp1 = ", cp1
		#	@info "cp2 = ", cp2
		#	@info "out = ", out
		# end
	end

	return take_log ? out : exp(out)
end

function rand(s::Stepdown; n = 1, take_log = false)
	u = rand(Float64, n)
	return map(p -> quantile(s, p, take_log = take_log), u)
end

function pdf(s::Stepdown, log_x; normalize = true, take_log = false)
	# Get the idx such that h_vals[idx] <= log(x) < h_vals[idx+1]
	idx = find_interval(log_x, s.log_x_vals)
	out = -Inf
	if idx <= s.N + 2
		out = s.log_h_vals[idx]
	end
	out = s.log_h_vals[idx]
	out = normalize ? out - s.norm_const : out
	return take_log ? out : exp(out)
end

# function pdf(s::Stepdown, x; normalize = true)
#	exp(logpdf(s, x; normalize))
# end

function cdf(s::Stepdown, log_x)
	idx = find_interval(log_x, s.log_x_vals)

	if idx == 0
		out = 0
	elseif idx >= s.N + 1
		out = 1
	else
		# Get the idx such that h_vals[idx] <= log(x) < h_vals[idx+1]	
		j1 = idx
		j2 = idx + 1
		cp1 = s.cum_probs[j1]
		cp2 = s.cum_probs[j2]
		log_x1 = s.log_x_vals[j1]
		log_x2 = s.log_x_vals[j2]

		# Do the following computation, but be careful to keep x values on
		# the log-scale, since they may be extremely small.
		# out = cp1 + (x - x1) / (x2 - x1) * (cp2 - cp1)
		log_ratio = log_sub(log_x, log_x1) - log_sub(log_x2, log_x1) + log(cp2 - cp1)
		out = exp(log_add(log(cp1), log_ratio))
	end

	return out
end

# Update step function with new log_u and log_p values. We could consider using
# heaps for the internal data structures to avoid re-sorting.
function add(s::Stepdown, log_u)
	log_h = log_p(s.w, s.g, log_u)
		
	# Add new x and h value, and keep track of the order in which it was added
	log_x_vals = cat(s.log_x_vals, log_u, dims = (1,1))
	log_h_vals = cat(s.log_h_vals, log_h, dims = (1,1))
	knot_order = cat(s.knot_order, s.N + 3, dims = (1,1))

	idx = sortperm(log_x_vals)
	s.log_x_vals = log_x_vals[idx]
	s.log_h_vals = log_h_vals[idx]
	s.knot_order = knot_order[idx]

	# Update cumulative probabilities
	s.N += 1
	(areas, cum_probs, norm_const) = steps2probs(s.log_x_vals, s.log_h_vals)

	return s
end

# Compute normalizing constant and cumulative probabilities based on
# log_x_vals and log_h_vals.
function steps2probs(log_x_vals, log_h_vals)
	if length(log_x_vals) != length(log_h_vals)
		error("log_x_vals and log_h_vals must be the same length")
	end

	N = length(log_x_vals) - 2
	areas = Vector{Float64}(undef, N+1)
	cum_probs = Vector{Float64}(undef, N+1)
	norm_const = 0

	if false
		# This is a more readable version of the calculation that may be less precise
		widths = diff(exp.(log_x_vals))
		heights = exp.(first(log_h_vals, N+1))
		areas = widths * heights
		norm_const = sum(areas)
		cum_probs = cumsum(areas) / norm_const
	else
		# This version tries to be more precise and do more work on the log-scale
		log_areas = zeros(N+1)
		log_cum_areas = zeros(N+1)

		log_areas[1] = log_h_vals[1] + log_x_vals[2]
		log_cum_areas[1] = log_areas[1]
		for l = 2:(N+1)
			log_areas[l] = log_h_vals[l] + log_sub(log_x_vals[l+1], log_x_vals[l])
			arg1 = max(log_cum_areas[l-1], log_areas[l])
			arg2 = min(log_cum_areas[l-1], log_areas[l])
			log_cum_areas[l] = log_add(arg1, arg2)
		end

		log_normconst = log_cum_areas[N+1]
		log_cum_probs = log_cum_areas .- log_normconst
		norm_const = exp.(log_normconst)
		cum_probs = exp.(log_cum_probs)
	end

	if any(map(isnan, cum_probs))
		error("NA values found in cum_probs. Approximation may have failed")
	end

	return areas, cum_probs, norm_const
end

