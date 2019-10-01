#include <fastenvelope/Logger.hpp>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/details/registry.h>
#include <spdlog/details/thread_pool.h>

#include <memory>
#include <mutex>
#include <iostream>

namespace fastEnvelope
{
std::shared_ptr<spdlog::async_logger> Logger::logger_;

// See https://github.com/gabime/spdlog#asynchronous-logger-with-multi-sinks
void Logger::init(std::vector<spdlog::sink_ptr> &sinks)
{
	auto l = spdlog::get("fastEnvelope");
	bool fastEnvelope = l != nullptr;
	if (fastEnvelope)
		spdlog::drop("fastEnvelope");

	spdlog::init_thread_pool(8192, 1);
	Logger::logger_ =
		std::make_shared<spdlog::async_logger>(
			"fastEnvelope",
			sinks.begin(), sinks.end(),
			spdlog::thread_pool(), spdlog::async_overflow_policy::block);
	spdlog::register_logger(logger_);

	if (fastEnvelope)
		logger().warn("Removed another fastEnvelope logger");
}

void Logger::init(bool use_cout, const std::string &filename, bool truncate)
{
	std::vector<spdlog::sink_ptr> sinks;
	if (use_cout)
	{
		sinks.emplace_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
	}
	if (!filename.empty())
	{
		sinks.emplace_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>(filename, truncate));
	}

	init(sinks);
}

void Logger::init(std::ostream &os)
{
	std::vector<spdlog::sink_ptr> sinks;
	sinks.emplace_back(std::make_shared<spdlog::sinks::ostream_sink_mt>(os, false));

	init(sinks);
}

} // namespace fastEnvelope
