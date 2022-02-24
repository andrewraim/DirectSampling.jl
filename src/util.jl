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

