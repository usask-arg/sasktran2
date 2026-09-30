#include "../../successive_orders/first_order_scratch.h"

#include <sasktran2/test_helper.h>

#include <condition_variable>
#include <future>
#include <mutex>
#include <stdexcept>

namespace {
    using sasktran2::successive_orders::FirstOrderScratchLease;
}

TEST_CASE("First-order scratch leases retain capacity and expose active spans",
          "[successive_orders][first_order_scratch]") {
    std::uint64_t allocation_id;
    const double* solar_address;
    const double* table_address;
    const double* endpoint_address;
    std::size_t payload_bytes;
    {
        FirstOrderScratchLease scratch(3, 32, 16);
        scratch.prepare_endpoints(64);
        REQUIRE(scratch.shared());
        allocation_id = scratch.allocation_id();
        auto solar = scratch.solar(0);
        auto table = scratch.table();
        auto extinction = scratch.endpoint_extinction();
        scratch.solar(1).setConstant(7.0);
        scratch.solar(2).setConstant(-3.0);
        solar.setConstant(2.0);
        table.setConstant(5.0);
        extinction.setConstant(-11.0);
        solar_address = solar.data();
        table_address = table.data();
        endpoint_address = extinction.data();
        payload_bytes = scratch.storage_bytes();
        REQUIRE(payload_bytes >= (3 * 32 + 16 + 2 * 64) * sizeof(double));
        REQUIRE_THROWS_AS(scratch.solar(-1), std::out_of_range);
        REQUIRE_THROWS_AS(scratch.solar(3), std::out_of_range);
    }
    {
        FirstOrderScratchLease scratch(1, 5, 3);
        scratch.prepare_endpoints(4);
        REQUIRE(scratch.allocation_id() == allocation_id);
        REQUIRE(scratch.solar(0).data() == solar_address);
        REQUIRE(scratch.table().data() == table_address);
        REQUIRE(scratch.endpoint_extinction().data() == endpoint_address);
        REQUIRE(scratch.solar(0).size() == 5);
        REQUIRE(scratch.table().size() == 3);
        REQUIRE(scratch.endpoint_extinction().size() == 4);
        REQUIRE(scratch.endpoint_albedo().size() == 4);
        REQUIRE(scratch.storage_bytes() == payload_bytes);
        scratch.solar(0).setZero();
        scratch.table().setZero();
        scratch.endpoint_extinction().setZero();
    }
    {
        FirstOrderScratchLease scratch(3, 32, 16);
        scratch.prepare_endpoints(64);
        REQUIRE(scratch.solar(0).head(5).isZero(0.0));
        REQUIRE(scratch.solar(0).tail(27).isConstant(2.0));
        REQUIRE(scratch.solar(1).isConstant(7.0));
        REQUIRE(scratch.solar(2).isConstant(-3.0));
        REQUIRE(scratch.table().tail(13).isConstant(5.0));
        REQUIRE(scratch.endpoint_extinction().tail(60).isConstant(-11.0));
    }
}

TEST_CASE("Reentrant first-order scratch leases cannot invalidate outer views",
          "[successive_orders][first_order_scratch]") {
    FirstOrderScratchLease outer(1, 9, 7);
    outer.prepare_endpoints(5);
    auto solar = outer.solar(0);
    auto table = outer.table();
    auto extinction = outer.endpoint_extinction();
    solar.setConstant(0.125);
    table.setConstant(-0.375);
    extinction.setConstant(0.625);
    const auto original_id = outer.allocation_id();
    const auto original_bytes = outer.storage_bytes();
    {
        FirstOrderScratchLease nested(4, 101, 103);
        nested.prepare_endpoints(107);
        REQUIRE_FALSE(nested.shared());
        REQUIRE(nested.allocation_id() != original_id);
        REQUIRE(nested.solar(0).data() != solar.data());
        nested.solar(0).setConstant(99.0);
        nested.table().setConstant(101.0);
        nested.endpoint_extinction().setConstant(103.0);
        REQUIRE_THROWS_AS(FirstOrderScratchLease(-1, 0, 0),
                          std::invalid_argument);
    }
    REQUIRE(outer.shared());
    REQUIRE(outer.allocation_id() == original_id);
    REQUIRE(outer.storage_bytes() == original_bytes);
    REQUIRE(solar.isConstant(0.125));
    REQUIRE(table.isConstant(-0.375));
    REQUIRE(extinction.isConstant(0.625));
    FirstOrderScratchLease next_nested(0, 0, 0);
    REQUIRE_FALSE(next_nested.shared());
    REQUIRE(next_nested.table().size() == 0);
}

TEST_CASE("First-order scratch leases release occupancy during unwinding",
          "[successive_orders][first_order_scratch]") {
    std::uint64_t original_id = 0;
    try {
        FirstOrderScratchLease scratch(1, 2, 3);
        original_id = scratch.allocation_id();
        REQUIRE(scratch.shared());
        throw std::runtime_error("Product interrupted");
    } catch (const std::runtime_error&) {
    }
    FirstOrderScratchLease scratch(1, 1, 1);
    REQUIRE(scratch.shared());
    REQUIRE(scratch.allocation_id() == original_id);
}

TEST_CASE("Concurrent calling threads own distinct first-order scratch arenas",
          "[successive_orders][first_order_scratch]") {
    std::mutex mutex;
    std::condition_variable ready;
    int arrivals = 0;
    const auto product = [&](double value) {
        FirstOrderScratchLease scratch(2, 17, 11);
        scratch.prepare_endpoints(13);
        scratch.solar(0).setConstant(value);
        scratch.solar(1).setConstant(value + 1.0);
        scratch.table().setConstant(value + 2.0);
        scratch.endpoint_extinction().setConstant(value + 3.0);
        {
            std::unique_lock<std::mutex> lock(mutex);
            ++arrivals;
            ready.notify_all();
            ready.wait(lock, [&]() { return arrivals == 2; });
        }
        return std::make_pair(
            scratch.allocation_id(),
            scratch.shared() && scratch.solar(0).isConstant(value) &&
                scratch.solar(1).isConstant(value + 1.0) &&
                scratch.table().isConstant(value + 2.0) &&
                scratch.endpoint_extinction().isConstant(value + 3.0));
    };
    auto first = std::async(std::launch::async, product, 3.0);
    auto second = std::async(std::launch::async, product, -7.0);
    const auto first_result = first.get();
    const auto second_result = second.get();
    REQUIRE(first_result.first != second_result.first);
    REQUIRE(first_result.second);
    REQUIRE(second_result.second);
}
