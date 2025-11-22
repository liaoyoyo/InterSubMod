#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <mutex>
#include <sstream>
#include <vector>

namespace InterSubMod {
namespace Utils {

enum class LogLevel {
    L_DEBUG,
    L_INFO,
    L_WARNING,
    L_ERROR
};

/**
 * @brief specific logger class that supports console and file output with different log levels.
 */
class Logger {
public:
    static Logger& instance();

    void set_log_level(LogLevel level);
    void set_log_file(const std::string& filename);
    
    void log(LogLevel level, const std::string& message);
    
    // Template for variadic arguments formatting if needed, 
    // but keeping it simple with string for now or use a helper.
    
    template<typename T>
    Logger& operator<<(const T& msg) {
        std::lock_guard<std::mutex> lock(mutex_);
        buffer_ << msg;
        return *this;
    }

    // Flush buffer with a specific level
    void end_log(LogLevel level);

    // Helper methods for quick logging
    static void debug(const std::string& msg);
    static void info(const std::string& msg);
    static void warning(const std::string& msg);
    static void error(const std::string& msg);

private:
    Logger() = default;
    ~Logger();
    
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    LogLevel current_level_ = LogLevel::L_INFO;
    std::ofstream log_file_;
    std::mutex mutex_;
    std::stringstream buffer_;
    
    std::string level_to_string(LogLevel level);
};

} // namespace Utils
} // namespace InterSubMod

