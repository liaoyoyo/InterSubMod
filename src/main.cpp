#include <iostream>
#include "core/Config.hpp"
#include "utils/ArgParser.hpp"

// Placeholder for future core logic includes
// #include "core/SomaticSnv.hpp"

int main(int argc, char** argv) {
    InterSubMod::Config config;
    
    if (!InterSubMod::Utils::ArgParser::parse(argc, argv, config)) {
        return 1; // Parse failed or help printed
    }

    std::cout << "Validating configuration..." << std::endl;
    if (!config.validate()) {
        std::cerr << "Configuration validation failed." << std::endl;
        return 1;
    }

    config.print();

    std::cout << "Configuration valid. Starting analysis (Not implemented yet)..." << std::endl;
    
    // Future main loop:
    // 1. Load VCF -> SomaticSnvTable
    // 2. Generate Regions
    // 3. For each region:
    //    Fetch Reads -> Parse MM -> Build Matrix -> Calculate Distance -> Cluster -> Assoc
    
    return 0;
}

