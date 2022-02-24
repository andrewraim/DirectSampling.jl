module DirectSampling

using DataStructures
using Dates
using Logging
using LoggingExtras
using Printf
using Random

import Base.rand

include("bisection.jl")
include("direct-sampler.jl")
include("find-interval.jl")
include("util.jl")
include("BaseDistribution.jl")
include("WeightFunction.jl")
include("Interval.jl")
include("Stepdown.jl")

export
	BaseDistribution
	Stepdown
	WeightFunction

	bisection
	direct_sampler
	find_interval
	log_c
	log_eval
	roots
end
