#include "utils/Logger.hpp"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <filesystem>

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
    
    // Ensure directory exists
    std::filesystem::path p(filename);
    if (p.has_parent_path()) {
        std::filesystem::create_directories(p.parent_path());
    }

    log_file_.open(filename, std::ios::app);
}

std::string Logger::level_to_string(LogLevel level) {
    switch (level) {
        case LogLevel::LOG_DEBUG: return "DEBUG";
        case LogLevel::LOG_INFO:  return "INFO ";
        case LogLevel::LOG_WARN:  return "WARN ";
        case LogLevel::LOG_ERROR: return "ERROR";
        default: return "UNK  ";
    }
}

std::string Logger::get_color_code(LogLevel level) {
    // ANSI color codes
    switch (level) {
        case LogLevel::LOG_DEBUG: return "\033[36m"; // Cyan
        case LogLevel::LOG_INFO:  return "\033[32m"; // Green
        case LogLevel::LOG_WARN:  return "\033[33m"; // Yellow
        case LogLevel::LOG_ERROR: return "\033[31m"; // Red
        default: return "";
    }
}

std::string Logger::reset_color_code() {
    return "\033[0m";
}

void Logger::log(LogLevel level, const std::string& message, const char* file, int line) {
    // Check level before locking
    // Since LogLevel (Types.hpp) is verbosity based (DEBUG=3 > ERROR=0),
    // we print if level <= current_level_.
    // Wait, typical comparison:
    // If msg is DEBUG(3) and current is INFO(2). 3 <= 2 False. Don't print.
    // If msg is ERROR(0) and current is INFO(2). 0 <= 2 True. Print.
    // So:
    if (static_cast<int>(level) > static_cast<int>(current_level_)) {
        return;
    }

    std::lock_guard<std::mutex> lock(mutex_);

    // Get current time
    auto now = std::chrono::system_clock::now();
    auto now_time = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

    std::stringstream ss;
    struct tm* time_info = std::localtime(&now_time);
    
    // Format: [Time][Thread][Level] Message (File:Line)
    ss << "[" << std::put_time(time_info, "%Y-%m-%d %H:%M:%S") << "." << std::setfill('0') << std::setw(3) << ms.count() << "]";
    
#ifdef _OPENMP
    ss << "[T" << omp_get_thread_num() << "]";
#endif

    ss << "[" << level_to_string(level) << "] " << message;

    if (file && (level == LogLevel::LOG_DEBUG || level == LogLevel::LOG_ERROR)) {
        std::filesystem::path p(file);
        ss << " (" << p.filename().string() << ":" << line << ")";
    }
    
    ss << std::endl;

    // Output to console with colors
    std::cout << get_color_code(level) << ss.str() << reset_color_code() << std::flush;

    // Output to file (no colors)
    if (log_file_.is_open()) {
        log_file_ << ss.str() << std::flush;
    }
}

void Logger::debug(const std::string& msg, const char* file, int line) {
    instance().log(LogLevel::LOG_DEBUG, msg, file, line);
}

void Logger::info(const std::string& msg, const char* file, int line) {
    instance().log(LogLevel::LOG_INFO, msg, file, line);
}

void Logger::warning(const std::string& msg, const char* file, int line) {
    instance().log(LogLevel::LOG_WARN, msg, file, line);
}

void Logger::error(const std::string& msg, const char* file, int line) {
    instance().log(LogLevel::LOG_ERROR, msg, file, line);
}

// ScopedLogger Implementation
ScopedLogger::ScopedLogger(const std::string& action_name, LogLevel level) 
    : action_name_(action_name), level_(level), start_time_(std::chrono::steady_clock::now()) {
    Logger::instance().log(level_, "START: " + action_name_);
}

ScopedLogger::~ScopedLogger() {
    auto end_time = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time_).count();
    Logger::instance().log(level_, "DONE : " + action_name_ + " (" + std::to_string(duration) + " ms)");
}

}  // namespace Utils
}  // namespace InterSubMod
