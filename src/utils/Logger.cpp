#include "utils/Logger.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>

namespace InterSubMod {
namespace Utils {

Logger& Logger::instance() {
    static Logger instance;
    return instance;
}

Logger::~Logger() {
    if (log_file_.is_open()) {
        log_file_.close();
    }
}

void Logger::set_log_level(LogLevel level) {
    std::lock_guard<std::mutex> lock(mutex_);
    current_level_ = level;
}

void Logger::set_log_file(const std::string& filename) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (log_file_.is_open()) {
        log_file_.close();
    }
    log_file_.open(filename, std::ios::app);
}

std::string Logger::level_to_string(LogLevel level) {
    switch (level) {
        case LogLevel::L_DEBUG:   return "DEBUG";
        case LogLevel::L_INFO:    return "INFO";
        case LogLevel::L_WARNING: return "WARNING";
        case LogLevel::L_ERROR:   return "ERROR";
        default:                  return "UNKNOWN";
    }
}

void Logger::log(LogLevel level, const std::string& message) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (level < current_level_) {
        return;
    }

    // Get current time
    auto now = std::chrono::system_clock::now();
    auto now_time = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

    std::stringstream ss;
    // Note: std::put_time is not thread-safe on some old platforms but generally ok in modern C++.
    // Using localtime_r or localtime_s is better but std::put_time takes tm*.
    struct tm* time_info = std::localtime(&now_time);
    
    ss << "[" << std::put_time(time_info, "%Y-%m-%d %H:%M:%S") 
       << "." << std::setfill('0') << std::setw(3) << ms.count() << "] "
       << "[" << level_to_string(level) << "] "
       << message << std::endl;

    // Output to console
    std::cout << ss.str();
    
    // Output to file if open
    if (log_file_.is_open()) {
        log_file_ << ss.str();
        log_file_.flush();
    }
}

void Logger::end_log(LogLevel level) {
    std::lock_guard<std::mutex> lock(mutex_);
    std::string msg = buffer_.str();
    buffer_.str(""); // Clear buffer
    buffer_.clear();
    
    // Note: calling log() from here is tricky if we hold mutex, but log() also locks.
    // To avoid deadlock, we shouldn't call log() while holding mutex.
    // But we just cleared buffer. 
    // The correct way would be:
}

// Static helpers
void Logger::debug(const std::string& msg) {
    instance().log(LogLevel::L_DEBUG, msg);
}

void Logger::info(const std::string& msg) {
    instance().log(LogLevel::L_INFO, msg);
}

void Logger::warning(const std::string& msg) {
    instance().log(LogLevel::L_WARNING, msg);
}

void Logger::error(const std::string& msg) {
    instance().log(LogLevel::L_ERROR, msg);
}

} // namespace Utils
} // namespace InterSubMod
