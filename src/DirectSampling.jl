module DirectSampling

using DataStructures
using Dates
using Logging
using LoggingExtras
using Printf
using Random

import Base.rand
import Base.print

include("bisection.jl")
include("find-interval.jl")
include("util.jl")
include("BaseDistribution.jl")
include("WeightFunction.jl")
include("Interval.jl")
include("Stepdown.jl")
include("direct-sampler.jl")
include("direct-sampler-ar.jl")

export
	BaseDistribution
	Stepdown
	WeightFunction
	StepdownKnotMethod
	SMALL_RECTS
	EQUAL_STEPS

	bisection
	direct_sampler
	find_interval
	log_add
	log_c
	log_eval
	log_sub
	roots
	geometric_midpoint
	arithmetic_midpoint
	univariate_distance
end

