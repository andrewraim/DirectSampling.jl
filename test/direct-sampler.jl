using Plots
import DirectSampling: direct_sampler, direct_sampler_ar, SMALL_RECTS

include("LognormalWeightFunction.jl")
include("LaplaceBaseDistribution.jl")

z = -10
mu = 3.5
sigma2 = 8.5
lambda = 0.4

w = LognormalWeightFunction(z, mu, sigma2)
g = LaplaceBaseDistribution(lambda)

out = direct_sampler(w, g, n = 5000, N = 100)
histogram(out, normalize = :probability, bins = 30)

(x, rejections) = direct_sampler_ar(w, g, n = 5000, N = 10)
histogram(x, normalize = :probability, bins = 30)

