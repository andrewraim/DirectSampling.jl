using Plots

include("direct-sampler.jl")

z = -10
mu = 3.5
sigma2 = 8.5
rho = 0.01

w_obj = LognormalWeightFunction(z, mu, sigma2)
g_obj = DGeomBaseDistribution(rho)
out = direct_sampler(500000, w_obj, g_obj, tol = 1e-10, N = 100, fill_method = "small_jumps");

histogram(out, normalize = :probability)
