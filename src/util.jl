# TBD: Maybe dealing with the logging here is overstepping our boundaries

function elapsed_seconds(then, now)
	el = now - then
	el.value / 1000
end

# Define a logger that includes timestamps
# From https://github.com/oxinabox/LoggingExtras.jl
const date_format = "yyyy-mm-dd HH:MM:SS"

timestamp_logger(logger) = TransformerLogger(logger) do log
	merge(log, (; message = "$(Dates.format(now(), date_format)) $(log.message)"))
end

ConsoleLogger(stdout, Logging.Debug) |> timestamp_logger |> global_logger

# The following function helps us to use the property:
#  log(x + y) = log(x) + log(1 + y/x)
#             = log(x) + log1p(exp(log(y) - log(x))),
#  with log(x) and log(y) given as inputs, i.e. on the log-scale. When x and y
#  are of very different magnitudes, this is more stable when x is taken to be
#  the larger of the inputs. The most extreme case is when one of the inputs
#  might be -inf; in this case that input should be the second one.
#
#  https://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation
function log_add(logx, logy)
	if logx < 0 && logy < 0 && isinf(logx) && isinf(logy)
		# Is it possible to handle this case more naturally?
		# this is just computing log(0 + 0) = -Inf
		return(-Inf)
	end
	logx + log1p(exp(logy - logx))
end

# Same as log_add, but for subtraction
#  log(x - y) = log(x) + log(1 - y/x)
#             = log(x) + log1p(-exp(log(y) - log(x)))
function log_sub(logx, logy)
	if logx < 0 && logy < 0 && isinf(logx) && isinf(logy)
		# Is it possible to handle this case more naturally?
		# this is just computing log(0 - 0) = -Inf
		return(-Inf)
	end
	logx + log1p(-exp(logy - logx))
end

function geometric_midpoint(log_x, log_y, take_log = false)
	out = 1/2 * (log_x + log_y)
	take_log ? out : exp(out)
end

function arithmetic_midpoint(log_x, log_y, take_log = false)
	out = log(1/2) + log_y + log1p(exp(log_x - log_y))
	take_log ? out : exp(out)
end

function univariate_distance(log_x, log_y, take_log = false)
	out = logsub(log_y, log_x)
	take_log ? out : exp(out)
end

