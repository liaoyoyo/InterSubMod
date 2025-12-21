#pragma once

#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>

#include "core/Types.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace InterSubMod {
namespace Utils {

/**
 * @brief Singleton Logger class for consistent debug output.
 */
class Logger {
public:
    static Logger& instance();

    void set_log_level(LogLevel level);
    void set_log_file(const std::string& filename);

    // Core logging function
    void log(LogLevel level, const std::string& message, const char* file = nullptr, int line = -1);

    // Static helpers for cleaner syntax
    static void debug(const std::string& msg, const char* file = nullptr, int line = -1);
    static void info(const std::string& msg, const char* file = nullptr, int line = -1);
    static void warning(const std::string& msg, const char* file = nullptr, int line = -1);
    static void error(const std::string& msg, const char* file = nullptr, int line = -1);

private:
    Logger() = default;
    ~Logger();

    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    LogLevel current_level_ = LogLevel::LOG_INFO;
    std::ofstream log_file_;
    std::mutex mutex_;

    std::string level_to_string(LogLevel level);
    std::string get_color_code(LogLevel level);
    std::string reset_color_code();
};

/**
 * @brief RAII helper to log start and end of a scope/action.
 */
class ScopedLogger {
public:
    ScopedLogger(const std::string& action_name, LogLevel level = LogLevel::LOG_INFO);
    ~ScopedLogger();

private:
    std::string action_name_;
    LogLevel level_;
    std::chrono::steady_clock::time_point start_time_;
};

}  // namespace Utils
}  // namespace InterSubMod

// Macros to automatically capture file and line number
#define LOG_DEBUG(msg) InterSubMod::Utils::Logger::debug(msg, __FILE__, __LINE__)
#define LOG_INFO(msg) InterSubMod::Utils::Logger::info(msg, __FILE__, __LINE__)
#define LOG_WARNING(msg) InterSubMod::Utils::Logger::warning(msg, __FILE__, __LINE__)
#define LOG_ERROR(msg) InterSubMod::Utils::Logger::error(msg, __FILE__, __LINE__)
