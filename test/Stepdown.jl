using Plots
using StatsPlots
import DirectSampling: Stepdown, geometric_midpoint, StepdownKnotMethod,
	SMALL_RECTS, pdf, cdf, quantile, rand

include("LognormalWeightFunction.jl")
include("LaplaceBaseDistribution.jl")

z = 2
mu = 0
sigma2 = 1
lambda = 0.4

log_midpoint(log_x, log_y) = geometric_midpoint(log_x, log_y; take_log = true)

w = LognormalWeightFunction(z, mu, sigma2)
g = LaplaceBaseDistribution(lambda)
step = Stepdown(w, g; tol = 1e-10, N = 10, knot_method = SMALL_RECTS,
	priority_weight = 1/2, log_midpoint = log_midpoint)

x_seq = map(exp, step.log_x_vals)
y_seq = map(exp, step.log_h_vals)
plot(x_seq, y_seq, seriestype = :scatter)

x_seq = 0:0.01:1
y_seq = map(x -> pdf(step, log(x), normalize = true), x_seq)
plot(x_seq, y_seq, seriestype = :scatter)

x_seq = 0:0.01:1
y_seq = map(x -> cdf(step, log(x)), x_seq)
plot(x_seq, y_seq, seriestype = :scatter)

y_seq = 0:0.01:1
x_seq = map(y -> quantile(step, y), y_seq)
plot(x_seq, y_seq, seriestype = :scatter)

out = rand(step, n = 50000, take_log = false)
histogram(out, normalize = :probability)

