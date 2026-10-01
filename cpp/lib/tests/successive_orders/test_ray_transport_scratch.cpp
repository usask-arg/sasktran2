#include "../../successive_orders/ray_transport_scratch.h"

#include <sasktran2/test_helper.h>

#include <array>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <future>
#include <limits>
#include <mutex>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
    using namespace sasktran2::successive_orders;

    struct LeasedRayFixture {
        LeasedRayFixture(int ray_count, int layer_count)
            : rays(ray_count),
              atmosphere(
                  sasktran2::atmosphere::AtmosphereGridStorageFull<1>(2, 3, 0),
                  sasktran2::atmosphere::Surface<1>(2), true),
              albedo((Eigen::MatrixXd(3, 2) << 0.1, 0.2, 0.3, 0.4, 0.7, 0.8)
                         .finished()) {
            atmosphere.storage().total_extinction << 0.3, 0.1, 0.5, 0.4, 0.2,
                0.7;
            atmosphere.storage().ssa << 0.8, 0.35, 0.6, 0.55, 0.4, 0.75;
            atmosphere.surface().set_spatial_lambertian_albedo(albedo);
            std::vector<int> offsets(static_cast<std::size_t>(ray_count) + 1);
            std::vector<int> columns, row_columns;
            for (int row = 0; row < ray_count; ++row) {
                auto& ray = rays[row];
                for (int layer = 0; layer < layer_count; ++layer) {
                    const auto offset = static_cast<std::uint32_t>(2 * layer);
                    ray.layers.wide_values().push_back(
                        {offset, 2, offset, 2, offset, 2});
                    ray.atmosphere_weights.push_back({layer % 3, 0.25});
                    ray.atmosphere_weights.push_back({(layer + 1) % 3, 0.75});
                    ray.source_weights.wide_values().emplace_back(
                        (row + layer) % 4, 0.35);
                    ray.source_weights.wide_values().emplace_back(
                        (row + layer + 2) % 4, 0.65);
                    ray.optical_depth_indices.push_back(layer % 3);
                    ray.optical_depth_indices.push_back((layer + 1) % 3);
                    ray.optical_depth_weights.push_back(0.3 + 0.1 * layer);
                    ray.optical_depth_weights.push_back(0.8 - 0.05 * layer);
                }
                ray.ground_hit = layer_count != 0 && row % 2 == 0;
                if (ray.ground_hit) {
                    ray.ground_weights = {{0, 0.2}, {3, 0.8}};
                    ray.ground_horizontal_weights = {{0, 0.25}, {2, 0.75}};
                }
                offsets[row] = static_cast<int>(columns.size());
                ray.transport_value_offset = columns.size();
                compile_transport_row(ray, row_columns);
                columns.insert(columns.end(), row_columns.begin(),
                               row_columns.end());
                compact_ray_interpolation(ray);
            }
            offsets[ray_count] = static_cast<int>(columns.size());
            sparsity =
                TransportSparsity(4, std::move(offsets), std::move(columns));
        }

        RayTransportMap make_map() const { return {rays, sparsity}; }

        void update(int stage) {
            if (stage == 1) {
                atmosphere.surface().set_spatial_lambertian_albedo(0.91 *
                                                                   albedo);
            } else if (stage == 2) {
                atmosphere.storage().total_extinction *= 1.03;
                atmosphere.storage().ssa *= 0.97;
            } else if (stage == 3) {
                atmosphere.storage().total_extinction << 0.3, 0.1, 0.5, 0.4,
                    0.2, 0.7;
                atmosphere.storage().ssa << 0.8, 0.35, 0.6, 0.55, 0.4, 0.75;
                atmosphere.surface().set_spatial_lambertian_albedo(albedo);
            }
        }

        std::vector<RayInterpolation> rays;
        TransportSparsity sparsity;
        sasktran2::atmosphere::Atmosphere<1> atmosphere;
        Eigen::MatrixXd albedo;
    };

    template <typename Actual, typename Expected>
    void leased_ray_require_bits(const Actual& actual,
                                 const Expected& expected) {
        REQUIRE(actual.size() == expected.size());
        REQUIRE(actual.allFinite());
        REQUIRE(expected.allFinite());
        if (actual.size() != 0) {
            REQUIRE(std::memcmp(actual.data(), expected.data(),
                                actual.size() * sizeof(double)) == 0);
        }
    }

    std::array<Eigen::VectorXd, 5>
    leased_ray_snapshot(const RayTransportWorkspace& workspace) {
        return {workspace.optical_depth, workspace.albedo,
                workspace.transmission_before, workspace.source_fraction,
                workspace.factor_cotangent};
    }
} // namespace

TEST_CASE("Scalar LOS transport leases reuse owning scratch between calls",
          "[successive_orders][ray_transport_scratch]") {
    std::uint64_t id;
    const double *values, *optical_depth;
    std::size_t bytes;
    {
        ScalarRayTransportWorkspaceLease lease;
        REQUIRE(lease.shared());
        id = lease.allocation_id();
        lease.values().setConstant(19, 0.125);
        lease.workspace().resize(7);
        lease.workspace().optical_depth.setConstant(-0.375);
        lease.workspace().albedo.setConstant(0.625);
        values = lease.values().data();
        optical_depth = lease.workspace().optical_depth.data();
        bytes = lease.storage_bytes();
        REQUIRE(bytes == (19 + 5 * 7) * sizeof(double));
        REQUIRE(bytes == lease.value_bytes() + lease.workspace_bytes());
    }
    ScalarRayTransportWorkspaceLease lease;
    REQUIRE(lease.shared());
    REQUIRE(lease.allocation_id() == id);
    REQUIRE(lease.values().data() == values);
    REQUIRE(lease.workspace().optical_depth.data() == optical_depth);
    REQUIRE(lease.values().isConstant(0.125));
    REQUIRE(lease.workspace().optical_depth.isConstant(-0.375));
    REQUIRE(lease.workspace().albedo.isConstant(0.625));
    REQUIRE(lease.storage_bytes() == bytes);
}

TEST_CASE("LOS pooled values survive failed growth and remain reusable",
          "[successive_orders][ray_transport_scratch]") {
    const auto overflow_size = std::numeric_limits<Eigen::Index>::max();
    REQUIRE(static_cast<std::size_t>(overflow_size) >
            std::numeric_limits<std::size_t>::max() / sizeof(double));
    ScalarRayTransportWorkspaceLease lease;
    lease.prepare_values(7);
    lease.values().setConstant(0.625);
    const auto* values = lease.values().data();
    REQUIRE_THROWS_AS(lease.prepare_values(overflow_size), std::bad_alloc);
    REQUIRE(lease.values().size() == 7);
    REQUIRE(lease.values().data() == values);
    REQUIRE(lease.values().isConstant(0.625));
    lease.prepare_values(7);
    REQUIRE(lease.values().data() == values);
    lease.prepare_values(0);
    REQUIRE(lease.values().size() == 0);
    lease.prepare_values(3);
    lease.values().setConstant(-0.375);
    REQUIRE(lease.values().isConstant(-0.375));
}

TEST_CASE("Ray transport repairs every retained companion before VJP",
          "[successive_orders][ray_transport_scratch]") {
    LeasedRayFixture fixture(2, 3);
    const auto map = fixture.make_map();
    ScalarRayTransportWorkspaceLease lease;
    lease.workspace().resize(3);
    lease.workspace().optical_depth.setConstant(0.625);
    const auto* optical_depth = lease.workspace().optical_depth.data();
    const auto bytes = lease.workspace_bytes();
    REQUIRE_THROWS_AS(lease.workspace().resize(-1), std::invalid_argument);
    REQUIRE(lease.workspace_bytes() == bytes);
    REQUIRE(lease.workspace().optical_depth.data() == optical_depth);
    REQUIRE(lease.workspace().optical_depth.isConstant(0.625));
    lease.workspace().resize(3);
    REQUIRE(lease.workspace().optical_depth.data() == optical_depth);

    // The old single-vector guard did not repair these missing companions.
    lease.workspace().albedo.resize(0);
    lease.workspace().transmission_before.resize(0);
    lease.workspace().source_fraction.resize(0);
    lease.workspace().factor_cotangent.resize(0);
    const Eigen::VectorXd cotangent =
        Eigen::VectorXd::LinSpaced(map.sparsity().nonzeros(), -0.125, 0.375);
    Eigen::VectorXd gradient =
        Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
    Eigen::VectorXd independent_gradient = gradient;
    RayTransportWorkspace independent;
    map.accumulate_vjp(fixture.atmosphere, 0, cotangent, gradient,
                       lease.workspace());
    map.accumulate_vjp(fixture.atmosphere, 0, cotangent, independent_gradient,
                       independent);
    REQUIRE(lease.workspace_bytes() == bytes);
    leased_ray_require_bits(gradient, independent_gradient);
    const auto actual_buffers = leased_ray_snapshot(lease.workspace());
    const auto expected_buffers = leased_ray_snapshot(independent);
    for (std::size_t index = 0; index < actual_buffers.size(); ++index) {
        leased_ray_require_bits(actual_buffers[index], expected_buffers[index]);
    }
}

TEST_CASE("Scalar LOS scratch reuse preserves complete native ray products",
          "[successive_orders][ray_transport_scratch][linearization]") {
    std::uint64_t id = 0;
    for (const auto shape :
         {std::make_pair(1, 1), std::make_pair(3, 4), std::make_pair(2, 0),
          std::make_pair(0, 0), std::make_pair(1, 2), std::make_pair(3, 4)}) {
        CAPTURE(shape.first, shape.second);
        LeasedRayFixture fixture(shape.first, shape.second);
        const auto map = fixture.make_map();
        TransportOperator transport(map.sparsity());
        const Eigen::VectorXd state = Eigen::VectorXd::LinSpaced(4, -0.4, 0.6);
        const Eigen::VectorXd state_tangent =
            Eigen::VectorXd::LinSpaced(4, 0.08, -0.11);
        const Eigen::VectorXd native_tangent = Eigen::VectorXd::LinSpaced(
            fixture.atmosphere.num_deriv(), -0.07, 0.09);
        const Eigen::VectorXd ray_cotangent =
            Eigen::VectorXd::LinSpaced(shape.first, -0.21, 0.31);
        RayTransportWorkspace independent_workspace;
        std::array<Eigen::VectorXd, 2> initial_value, initial_jvp, initial_vjp;
        for (int stage = 0; stage < 4; ++stage) {
            CAPTURE(stage);
            fixture.update(stage);
            for (int wavelength = 0; wavelength < 2; ++wavelength) {
                CAPTURE(wavelength);
                map.assemble_values(fixture.atmosphere, wavelength, transport);
                Eigen::VectorXd ray_value(shape.first);
                transport.apply_stokes<1>(state, ray_value);
                Eigen::VectorXd independent_values(map.sparsity().nonzeros());
                map.assemble_jvp(fixture.atmosphere, wavelength, native_tangent,
                                 independent_values);
                Eigen::VectorXd independent_jvp(shape.first);
                transport.apply_jvp_stokes<1>(
                    state, state_tangent, independent_values, independent_jvp);
                Eigen::VectorXd pooled_jvp(shape.first), saved_value_tangent;
                {
                    ScalarRayTransportWorkspaceLease lease;
                    REQUIRE(lease.shared());
                    if (id == 0) {
                        id = lease.allocation_id();
                    }
                    REQUIRE(lease.allocation_id() == id);
                    lease.values().setConstant(map.sparsity().nonzeros(),
                                               123.0);
                    map.assemble_jvp(fixture.atmosphere, wavelength,
                                     native_tangent, lease.values());
                    leased_ray_require_bits(lease.values(), independent_values);
                    transport.apply_jvp_stokes<1>(state, state_tangent,
                                                  lease.values(), pooled_jvp);
                    saved_value_tangent = lease.values();
                }
                leased_ray_require_bits(pooled_jvp, independent_jvp);

                Eigen::VectorXd independent_state_cotangent(4),
                    pooled_state_cotangent(4);
                const Eigen::VectorXd gradient_seed =
                    Eigen::VectorXd::LinSpaced(fixture.atmosphere.num_deriv(),
                                               -0.13, 0.17);
                Eigen::VectorXd independent_gradient = gradient_seed;
                Eigen::VectorXd pooled_gradient = gradient_seed;
                transport.apply_vjp_stokes<1>(state, ray_cotangent,
                                              independent_state_cotangent,
                                              independent_values);
                map.accumulate_vjp(fixture.atmosphere, wavelength,
                                   independent_values, independent_gradient,
                                   independent_workspace);
                {
                    ScalarRayTransportWorkspaceLease lease;
                    REQUIRE(lease.allocation_id() == id);
                    lease.values().setConstant(map.sparsity().nonzeros(),
                                               -321.0);
                    transport.apply_vjp_stokes<1>(state, ray_cotangent,
                                                  pooled_state_cotangent,
                                                  lease.values());
                    leased_ray_require_bits(lease.values(), independent_values);
                    map.accumulate_vjp(fixture.atmosphere, wavelength,
                                       lease.values(), pooled_gradient,
                                       lease.workspace());
                }
                leased_ray_require_bits(pooled_state_cotangent,
                                        independent_state_cotangent);
                leased_ray_require_bits(pooled_gradient, independent_gradient);
                // These owned products survive the lease's reuse for VJP.
                leased_ray_require_bits(pooled_jvp, independent_jvp);
                Eigen::VectorXd restored_value_tangent(
                    map.sparsity().nonzeros());
                map.assemble_jvp(fixture.atmosphere, wavelength, native_tangent,
                                 restored_value_tangent);
                leased_ray_require_bits(saved_value_tangent,
                                        restored_value_tangent);
                if (stage == 0) {
                    initial_value[wavelength] = ray_value;
                    initial_jvp[wavelength] = pooled_jvp;
                    initial_vjp[wavelength] = pooled_gradient;
                } else if (stage == 3) {
                    leased_ray_require_bits(ray_value,
                                            initial_value[wavelength]);
                    leased_ray_require_bits(pooled_jvp,
                                            initial_jvp[wavelength]);
                    leased_ray_require_bits(pooled_gradient,
                                            initial_vjp[wavelength]);
                }
            }
        }
    }
}

TEST_CASE("Nested LOS transport products cannot overwrite the outer lease",
          "[successive_orders][ray_transport_scratch][reentrant]") {
    LeasedRayFixture fixture(3, 11);
    const auto map = fixture.make_map();
    const Eigen::VectorXd native_tangent =
        Eigen::VectorXd::LinSpaced(fixture.atmosphere.num_deriv(), -0.07, 0.09);
    ScalarRayTransportWorkspaceLease outer;
    outer.values().setConstant(5, 0.125);
    outer.workspace().resize(2);
    outer.workspace().optical_depth.setConstant(0.3);
    outer.workspace().albedo.setConstant(0.4);
    outer.workspace().transmission_before.setConstant(0.5);
    outer.workspace().source_fraction.setConstant(0.6);
    outer.workspace().factor_cotangent.setConstant(0.7);
    const auto snapshot = leased_ray_snapshot(outer.workspace());
    const double* outer_values = outer.values().data();
    {
        ScalarRayTransportWorkspaceLease nested;
        REQUIRE_FALSE(nested.shared());
        REQUIRE(nested.allocation_id() != outer.allocation_id());
        nested.values().resize(map.sparsity().nonzeros());
        map.assemble_jvp(fixture.atmosphere, 1, native_tangent,
                         nested.values());
        Eigen::VectorXd gradient =
            Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
        map.accumulate_vjp(fixture.atmosphere, 1, nested.values(), gradient,
                           nested.workspace());
        REQUIRE(nested.workspace().optical_depth.size() == 11);
        REQUIRE(gradient.allFinite());
    }
    REQUIRE(outer.values().data() == outer_values);
    REQUIRE(outer.values().isConstant(0.125));
    const auto after = leased_ray_snapshot(outer.workspace());
    for (std::size_t index = 0; index < snapshot.size(); ++index) {
        leased_ray_require_bits(after[index], snapshot[index]);
    }
}

TEST_CASE("LOS transport scratch leases release after native failures",
          "[successive_orders][ray_transport_scratch]") {
    LeasedRayFixture fixture(2, 3);
    const auto map = fixture.make_map();
    std::uint64_t id = 0;
    try {
        ScalarRayTransportWorkspaceLease lease;
        id = lease.allocation_id();
        lease.values().resize(map.sparsity().nonzeros());
        const Eigen::VectorXd invalid_tangent = Eigen::VectorXd::Zero(1);
        map.assemble_jvp(fixture.atmosphere, 0, invalid_tangent,
                         lease.values());
        FAIL("Invalid native tangent must throw");
    } catch (const std::invalid_argument&) {
    }
    ScalarRayTransportWorkspaceLease lease;
    REQUIRE(lease.shared());
    REQUIRE(lease.allocation_id() == id);
    lease.values().resize(map.sparsity().nonzeros());
    const Eigen::VectorXd tangent =
        Eigen::VectorXd::Ones(fixture.atmosphere.num_deriv());
    map.assemble_jvp(fixture.atmosphere, 1, tangent, lease.values());
    REQUIRE(lease.values().allFinite());
}

TEST_CASE("Concurrent OS workers own distinct LOS transport scratch",
          "[successive_orders][ray_transport_scratch]") {
    std::mutex mutex;
    std::condition_variable ready;
    int arrivals = 0;
    const auto calculate = [&](int layers) {
        LeasedRayFixture fixture(2, layers);
        const auto map = fixture.make_map();
        const Eigen::VectorXd tangent =
            Eigen::VectorXd::Ones(fixture.atmosphere.num_deriv());
        ScalarRayTransportWorkspaceLease lease;
        lease.values().resize(map.sparsity().nonzeros());
        map.assemble_jvp(fixture.atmosphere, 1, tangent, lease.values());
        Eigen::VectorXd gradient =
            Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
        map.accumulate_vjp(fixture.atmosphere, 1, lease.values(), gradient,
                           lease.workspace());
        const Eigen::VectorXd saved_values = lease.values();
        const auto saved_workspace = leased_ray_snapshot(lease.workspace());
        {
            std::unique_lock<std::mutex> lock(mutex);
            ++arrivals;
            ready.notify_all();
            ready.wait(lock, [&]() { return arrivals == 2; });
        }
        bool unchanged = lease.shared() && gradient.allFinite() &&
                         std::memcmp(saved_values.data(), lease.values().data(),
                                     saved_values.size() * sizeof(double)) == 0;
        const auto after = leased_ray_snapshot(lease.workspace());
        for (std::size_t index = 0; index < after.size(); ++index) {
            unchanged &=
                after[index].size() == saved_workspace[index].size() &&
                std::memcmp(after[index].data(), saved_workspace[index].data(),
                            after[index].size() * sizeof(double)) == 0;
        }
        return std::make_pair(lease.allocation_id(), unchanged);
    };
    auto first = std::async(std::launch::async, calculate, 3);
    auto second = std::async(std::launch::async, calculate, 7);
    const auto first_result = first.get();
    const auto second_result = second.get();
    REQUIRE(first_result.first != second_result.first);
    REQUIRE(first_result.second);
    REQUIRE(second_result.second);
}
