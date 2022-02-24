# using Plots
# using Gadfly
# using DirectSampling
# import DirectSampling
# import DirectSampling.Stepdown
import DirectSampling: Stepdown

include("LognormalWeightFunction.jl")
include("DGeomBaseDistribution.jl")

z = 2
mu = 0
sigma2 = 1
lambda = 0.4

w = LognormalWeightFunction(z, mu, sigma2)
g = LaplaceBaseDistribution(lambda)
step = Stepdown(w, g; tol = 1e-10, N = 10, method = "small_rects")

x_seq = map(exp, step.log_x_vals)
y_seq = map(exp, step.log_h_vals)
plot(x_seq, y_seq, seriestype = :scatter)

x_seq = 0:0.01:1
y_seq = map(x -> pdf(step, x; normalize = false), x_seq)
plot(x_seq, y_seq, seriestype = :scatter)

x_seq = 0:0.01:1
y_seq = map(x -> cdf(step, x), x_seq)
plot(x_seq, y_seq, seriestype = :scatter)

y_seq = 0:0.01:1
x_seq = map(y -> quantile(step, y), y_seq)
plot(x_seq, y_seq, seriestype = :scatter)

# TBD: Something seems to be going horribly wrong here, even though quantile
# function appears to be correct...
out = rand(step, 50000)
histogram(out, normalize = :probability)

let p
	p = Gadfly.plot(x = out);
	push!(p, Geom.density);
	p
	# p |> SVG("scatter.svg", 4inch, 4inch)
end

u = rand(Float64, 50)
Qu = map(x -> quantile(step,x), u)
histogram(Qu, normalize = :probability)



