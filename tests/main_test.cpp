#include <gtest/gtest.h>

#include "utils/ResourceMonitor.hpp"

using namespace InterSubMod::Utils;

class ResourceListener : public testing::EmptyTestEventListener {
public:
    void OnTestStart(const testing::TestInfo& /*test_info*/) override {
        monitor_.reset();
    }

    void OnTestEnd(const testing::TestInfo& test_info) override {
        std::string test_name = std::string(test_info.test_case_name()) + "." + test_info.name();
        monitor_.print_stats(test_name);
    }

private:
    ResourceMonitor monitor_;
};

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);

    testing::TestEventListeners& listeners = testing::UnitTest::GetInstance()->listeners();
    listeners.Append(new ResourceListener);

    return RUN_ALL_TESTS();
}
