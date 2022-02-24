"""
Define abstract WeightFunction

Catch-all methods for logeval, roots, and log_c are given here, but these must
be implemented specifically for each weight function.

TBD: some weight functions need state. Does this design work?
"""
abstract type WeightFunction end

function eval(w::T, x) where {T<:WeightFunction}
	exp(log_eval(w, x))
end

function log_eval(w::T, x) where {T<:WeightFunction}
	error("log_eval: please implement for your WeightFunction")
end

function roots(w::T, log_a) where {T<:WeightFunction}
	error("roots: please implement for your WeightFunction")
end

function log_c(w::T) where {T<:WeightFunction}
	error("log_c: please implement for your WeightFunction")
end

