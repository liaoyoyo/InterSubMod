#pragma once

#include <chrono>
#include <iostream>
#include <string>
#include <iomanip>

#ifdef USE_JEMALLOC
#include <jemalloc/jemalloc.h>
#endif

namespace InterSubMod {
namespace Utils {

class ResourceMonitor {
public:
    ResourceMonitor() {
        reset();
    }

    void reset() {
        start_time_ = std::chrono::high_resolution_clock::now();
    }

    double get_elapsed_seconds() const {
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time_;
        return elapsed.count();
    }

    // Returns allocated memory in bytes
    size_t get_memory_usage() const {
        size_t allocated = 0;
#ifdef USE_JEMALLOC
        size_t sz = sizeof(size_t);
        // epoch needs to be advanced to get up-to-date stats
        uint64_t epoch = 1;
        mallctl("epoch", &epoch, &sz, &epoch, sizeof(epoch));

        if (mallctl("stats.allocated", &allocated, &sz, NULL, 0) != 0) {
            // Fail silently or log?
        }
#endif
        return allocated;
    }
    
    void print_stats(const std::string& label = "Execution") const {
        double time = get_elapsed_seconds();
        size_t mem = get_memory_usage();
        
        std::cout << "[" << label << "] ";
        std::cout << "Time: " << std::fixed << std::setprecision(4) << time << " s";
        
#ifdef USE_JEMALLOC
        std::cout << ", Memory: " << std::fixed << std::setprecision(2) << (mem / 1024.0 / 1024.0) << " MB";
#else
        std::cout << " (jemalloc not enabled)";
#endif
        std::cout << std::endl;
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;
};

} // namespace Utils
} // namespace InterSubMod

