#include <sasktran2/test_helper.h>
#include <sasktran2/wavelength_block.h>

#include <sasktran2.h>

TEST_CASE("Derivative-free wavelength block lane views", "[dual][batch]") {
    sasktran2::WavelengthBlockDual<3> block;
    block.resize(4, 0, true);

    for (int lane = 0; lane < block.block_capacity(); ++lane) {
        sasktran2::WavelengthBlockLaneDualView<3> lane_view(block, lane);
        REQUIRE(lane_view.derivative_size() == 0);
        REQUIRE(lane_view.deriv.size() == 0);
        lane_view.value.setConstant(static_cast<double>(lane + 1));
    }

    const auto& const_block = block;
    for (int lane = 0; lane < block.block_capacity(); ++lane) {
        const sasktran2::WavelengthBlockConstLaneDualView<3> lane_view(
            const_block, lane);
        REQUIRE(lane_view.deriv.size() == 0);
        REQUIRE(lane_view.value.isConstant(static_cast<double>(lane + 1)));
    }
}
