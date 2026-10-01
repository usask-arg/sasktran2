#include "../../successive_orders/first_order.h"

#include <sasktran2/test_helper.h>

#include <array>
#include <cfenv>
#include <cstring>
#include <memory>

#ifdef SKTRAN_RUST_SUPPORT

namespace {
    using namespace sasktran2::successive_orders;
    using Atmosphere = sasktran2::atmosphere::Atmosphere<1>;
    constexpr int wavelengths = 2;
    constexpr int moments = 8;

    class FirstOrderTestRounding {
      public:
        FirstOrderTestRounding() : m_original(std::fegetround()) {
            REQUIRE(m_original != -1);
            REQUIRE(std::fesetround(FE_TONEAREST) == 0);
        }
        ~FirstOrderTestRounding() { std::fesetround(m_original); }

      private:
        int m_original;
    };

    sasktran2::Config configuration(int threads = 1) {
        sasktran2::Config config;
        config.set_num_threads(threads);
        config.set_threading_model(
            sasktran2::Config::ThreadingModel::wavelength);
        config.set_num_do_streams(moments);
        config.set_num_singlescatter_moments(moments);
        config.set_num_hr_incoming(6);
        config.set_num_hr_outgoing(6);
        config.set_single_scatter_source(
            sasktran2::Config::SingleScatterSource::exact);
        config.set_multiple_scatter_source(
            sasktran2::Config::MultipleScatterSource::successive_orders);
        config.set_apply_delta_scaling(false);
        return config;
    }

    sasktran2::Geometry2D geometry() {
        Eigen::VectorXd altitudes(3);
        altitudes << 0.0, 1000.0, 3000.0;
        Eigen::VectorXd angles(3);
        angles << -0.02, 0.0, 0.02;
        return sasktran2::Geometry2D(0.6, 0.2, 6372000.0, std::move(altitudes),
                                     std::move(angles));
    }

    SourceGeometrySettings source_settings() {
        SourceGeometrySettings result;
        result.num_incoming = 6;
        result.num_outgoing = 6;
        result.num_sza = 2;
        result.num_threads = 1;
        result.altitude_grid_m = {500.0, 1750.0};
        return result;
    }

    Eigen::MatrixXd surface_albedo(const sasktran2::Geometry2D& geo) {
        Eigen::MatrixXd result(geo.num_horizontal_locations(), wavelengths);
        for (int horizontal = 0; horizontal < result.rows(); ++horizontal) {
            for (int wavelength = 0; wavelength < wavelengths; ++wavelength) {
                result(horizontal, wavelength) =
                    0.12 + 0.03 * horizontal + 0.02 * wavelength;
            }
        }
        return result;
    }

    void initialize_atmosphere(Atmosphere& atmosphere,
                               const sasktran2::Geometry2D& geo, int groups,
                               bool swapped_phase = false,
                               bool primal_degree_three = true) {
        auto& storage = atmosphere.storage();
        storage.leg_coeff.setZero();
        for (int wavelength = 0; wavelength < wavelengths; ++wavelength) {
            const bool uniform = wavelength == (swapped_phase ? 1 : 0);
            for (int location = 0; location < geo.size(); ++location) {
                storage.total_extinction(location, wavelength) =
                    (1.4e-5 + 1.1e-6 * location) * (0.9 + 0.2 * wavelength);
                storage.ssa(location, wavelength) =
                    0.87 - 0.005 * location - 0.01 * wavelength;
                storage.leg_coeff(0, location, wavelength) = 1.0;
                storage.leg_coeff(1, location, wavelength) =
                    0.08 + (uniform ? 0.0 : 0.004 * location);
                storage.leg_coeff(2, location, wavelength) =
                    0.35 + (uniform ? 0.0 : 0.007 * location);
                if (primal_degree_three) {
                    storage.leg_coeff(3, location, wavelength) =
                        0.025 + (uniform ? 0.0 : 0.0005 * location);
                }
            }
        }
        atmosphere.surface().set_spatial_lambertian_albedo(surface_albedo(geo));
        storage.resize_derivatives(groups);
        // A native volume VJP is requested even when phase groups are absent.
        auto& mapping = storage.get_derivative_mapping("native_volume");
        mapping.allocate_extinction_derivatives();
        mapping.allocate_ssa_derivatives();
        mapping.native_mapping().d_extinction->setOnes();
        mapping.native_mapping().d_ssa->setOnes();
        if (groups != 0) {
            mapping.allocate_legendre_derivatives();
            mapping.native_mapping().d_legendre->setZero();
            mapping.native_mapping().d_legendre->chip(1, 0).setConstant(0.015);
            // The separate adversarial case omits this degree from the primal.
            mapping.native_mapping().d_legendre->chip(3, 0).setConstant(0.03);
            mapping.native_mapping().scat_factor->setOnes();
            storage.finalize_scattering_derivatives(0);
        }
        storage.determine_maximum_order();
        REQUIRE(atmosphere.num_scattering_deriv_groups() == groups);
        REQUIRE(storage.max_order(0, 0) == (primal_degree_three ? 4 : 3));
        if (groups > 0) {
            REQUIRE(storage.d_max_order[0](0, 0) == 4);
        }
    }

    struct ProviderPair {
        ProviderPair()
            : config(configuration()), geo(geometry()), tracer(geo),
              captured_geometry(tracer, geo), original_geometry(tracer, geo),
              captured(geo, tracer), original(geo, tracer) {
            const sasktran2::viewinggeometry::InternalViewingGeometry viewing;
            captured.initialize_config(config);
            original.initialize_config(config);
            REQUIRE(captured.requests_endpoint_factors());
            REQUIRE(original.requests_endpoint_factors());
            captured_geometry.initialize(viewing, source_settings(), true);
            original_geometry.initialize(viewing, source_settings(), false);
            REQUIRE(!captured_geometry.incoming_endpoint_factors().empty());
            REQUIRE(original_geometry.incoming_endpoint_factors().empty());
            captured.initialize_geometry(captured_geometry);
            original.initialize_geometry(original_geometry);
            REQUIRE(captured.uses_compact_scalar_kernel());
            REQUIRE(original.uses_compact_scalar_kernel());
            REQUIRE(captured.can_release_incoming_geometry());
            REQUIRE(original.can_release_incoming_geometry());
            captured_geometry.release_incoming_traced_rays();
            original_geometry.release_incoming_traced_rays();
            REQUIRE(captured_geometry.incoming_endpoint_factors().empty());
            REQUIRE(captured_geometry.incoming_endpoint_factors().capacity() ==
                    0);
        }

        void initialize(Atmosphere& atmosphere, bool volume_changed = true) {
            captured.initialize_atmosphere(atmosphere, volume_changed);
            original.initialize_atmosphere(atmosphere, volume_changed);
        }

        sasktran2::Config config;
        sasktran2::Geometry2D geo;
        sasktran2::raytracing::RustRayTracer2D tracer;
        SourceGeometry1D captured_geometry, original_geometry;
        FirstOrderProvider<1> captured, original;
    };

    template <typename Actual, typename Expected>
    void require_bits(const Actual& actual, const Expected& expected) {
        REQUIRE(actual.size() == expected.size());
        REQUIRE(actual.allFinite());
        REQUIRE(expected.allFinite());
        if (actual.size() != 0) {
            REQUIRE(std::memcmp(actual.data(), expected.data(),
                                actual.size() * sizeof(double)) == 0);
        }
    }

    struct Products : std::array<Eigen::VectorXd, 13> {
        // Preserve each dot value and the residual, including the known
        // uncaptured-path discrepancy for derivative-only higher orders.
        std::array<double, 3> adjoint_dot{};
    };

    Products products(FirstOrderProvider<1>& provider,
                      const SourceGeometry1D& source,
                      const Atmosphere& atmosphere, int wavelength,
                      bool check_adjoint = true) {
        const int locations = atmosphere.storage().total_extinction.rows();
        Eigen::VectorXd tangent = Eigen::VectorXd::Zero(atmosphere.num_deriv());
        tangent.head(locations) =
            Eigen::VectorXd::LinSpaced(locations, -1.0e-7, 2.0e-7);
        tangent.segment(atmosphere.ssa_deriv_start_index(), locations) =
            Eigen::VectorXd::LinSpaced(locations, -0.02, 0.03);
        if (atmosphere.num_scattering_deriv_groups() > 0) {
            tangent.segment(atmosphere.scat_deriv_start_index(), locations) =
                Eigen::VectorXd::LinSpaced(locations, 0.04, -0.01);
        }
        tangent.tail(atmosphere.surface().num_deriv()).setConstant(0.025);
        const Eigen::VectorXd state =
            Eigen::VectorXd::LinSpaced(source.total_num_outgoing(), -0.1, 0.4);
        const Eigen::VectorXd cotangent =
            Eigen::VectorXd::LinSpaced(provider.size(), -0.2, 0.5);
        Products result;
        for (std::size_t index = 0; index < 8; ++index) {
            result[index].resize(provider.size());
        }
        for (std::size_t index = 8; index < 11; ++index) {
            result[index].setZero(atmosphere.num_deriv());
        }
        provider.prepare_wavelength(wavelength, 0);
        provider.calculate(wavelength, 0, result[0]);
        TransportOperator transport(source.transport_sparsity());
        provider.calculate_with_transport(wavelength, 0, transport, result[1]);
        result[2] = transport.values();
        provider.project_transport_state(state, result[11], result[12]);
        // Compact JVP preserves the caller's already prepared primal forcing.
        result[3] = result[0];
        provider.calculate_jvp(wavelength, 0, tangent, result[3], result[4]);
        require_bits(result[3], result[0]);
        result[5] = result[1];
        provider.calculate_jvp_with_transport(wavelength, 0, tangent,
                                              result[11], result[12], result[7],
                                              result[5], result[6]);
        require_bits(result[5], result[1]);
        provider.accumulate_vjp(wavelength, 0, cotangent, result[8]);
        provider.accumulate_vjp_with_transport(wavelength, 0, state, cotangent,
                                               result[9]);
        provider.accumulate_vjp_with_projected_transport(
            wavelength, 0, result[11], result[12], cotangent, result[10]);
        REQUIRE(result[8].head(locations).cwiseAbs().maxCoeff() > 0.0);
        result.adjoint_dot = {tangent.dot(result[8]), cotangent.dot(result[4]),
                              0.0};
        result.adjoint_dot[2] = result.adjoint_dot[0] - result.adjoint_dot[1];
        if (check_adjoint) {
            REQUIRE(result.adjoint_dot[0] ==
                    Catch::Approx(result.adjoint_dot[1])
                        .epsilon(1.0e-10)
                        .margin(1.0e-12));
        }
        return result;
    }

    void require_products_equal(const Products& actual,
                                const Products& expected) {
        for (std::size_t index = 0; index < actual.size(); ++index) {
            CAPTURE(index);
            require_bits(actual[index], expected[index]);
        }
        REQUIRE(std::memcmp(actual.adjoint_dot.data(),
                            expected.adjoint_dot.data(),
                            actual.adjoint_dot.size() * sizeof(double)) == 0);
    }
} // namespace

namespace {
    void require_updated_products(int groups, bool primal_degree_three) {
        const bool check_adjoint = groups == 0 || primal_degree_three;
        ProviderPair pair;
        auto atmosphere = std::make_unique<Atmosphere>(wavelengths, pair.geo,
                                                       pair.config, true);
        initialize_atmosphere(*atmosphere, pair.geo, groups, false,
                              primal_degree_three);
        pair.initialize(*atmosphere);
        REQUIRE(pair.captured.workspace_bytes() <
                pair.original.workspace_bytes());
        std::array<Products, wavelengths> initial;
        for (int stage = 0; stage < 5; ++stage) {
            CAPTURE(groups, stage);
            const auto volume_revision = atmosphere->volume_revision();
            if (stage == 1) {
                atmosphere->surface().set_spatial_lambertian_albedo(
                    1.2 * surface_albedo(pair.geo));
                atmosphere->mark_surface_changed();
                REQUIRE(atmosphere->volume_revision() == volume_revision);
            } else if (stage == 2) {
                initialize_atmosphere(*atmosphere, pair.geo, groups, true,
                                      primal_degree_three);
                atmosphere->storage().total_extinction *= 1.15;
                atmosphere->storage().ssa.array() -= 0.03;
                atmosphere->mark_changed();
                REQUIRE(atmosphere->volume_revision() > volume_revision);
            } else if (stage == 3) {
                initialize_atmosphere(*atmosphere, pair.geo, groups, false,
                                      primal_degree_three);
                atmosphere->mark_changed();
            } else if (stage == 4) {
                const auto instance = atmosphere->instance_id();
                atmosphere = std::make_unique<Atmosphere>(wavelengths, pair.geo,
                                                          pair.config, true);
                initialize_atmosphere(*atmosphere, pair.geo, groups, false,
                                      primal_degree_three);
                REQUIRE(atmosphere->instance_id() != instance);
            }
            pair.initialize(*atmosphere, stage != 1);
            std::array<Products, wavelengths> current;
            for (int wavelength = 0; wavelength < wavelengths; ++wavelength) {
                CAPTURE(wavelength);
                current[wavelength] =
                    products(pair.captured, pair.captured_geometry, *atmosphere,
                             wavelength, check_adjoint);
                const auto original =
                    products(pair.original, pair.original_geometry, *atmosphere,
                             wavelength, check_adjoint);
                require_products_equal(current[wavelength], original);
                if (!primal_degree_three && stage == 0 && wavelength == 0) {
                    // Existing VJP endpoint truncation follows primal
                    // max_order, while JVP follows derivative max_order.
                    // Preserve that inherited discrepancy rather than changing
                    // native results.
                    REQUIRE(std::abs(current[wavelength].adjoint_dot[2]) >
                            1.0e-12);
                }
                if (stage == 0) {
                    initial[wavelength] = current[wavelength];
                } else if (stage >= 3) {
                    require_products_equal(current[wavelength],
                                           initial[wavelength]);
                }
            }
            // Wavelength one evicts zero's physical caches on this single
            // worker.
            for (int wavelength = 0; wavelength < wavelengths; ++wavelength) {
                const auto repeated =
                    products(pair.captured, pair.captured_geometry, *atmosphere,
                             wavelength, check_adjoint);
                const auto original =
                    products(pair.original, pair.original_geometry, *atmosphere,
                             wavelength, check_adjoint);
                require_products_equal(repeated, current[wavelength]);
                require_products_equal(repeated, original);
            }
        }
    }
} // namespace

TEST_CASE("Factored first-order endpoints preserve complete native products "
          "after spectral evictions and atmosphere updates",
          "[successive_orders][first_order][endpoint_factors][linearization]") {
    FirstOrderTestRounding rounding;
    const int groups = GENERATE(0, 1);
    require_updated_products(groups, true);
}

TEST_CASE("Factored first-order endpoints preserve the uncaptured higher-order "
          "phase derivative products and dot residual",
          "[successive_orders][first_order][endpoint_factors][linearization]") {
    FirstOrderTestRounding rounding;
    require_updated_products(1, false);
}

TEST_CASE("First-order rounding changes materialize original endpoint "
          "coefficients before physical calculations",
          "[successive_orders][first_order][endpoint_factors][rounding]") {
    FirstOrderTestRounding rounding;
    const int mode = GENERATE(FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO);
    ProviderPair pair;
    Atmosphere atmosphere(wavelengths, pair.geo, pair.config, true);
    initialize_atmosphere(atmosphere, pair.geo, 1);
    pair.initialize(atmosphere);
    REQUIRE(pair.captured.workspace_bytes() < pair.original.workspace_bytes());
    REQUIRE(std::fesetround(mode) == 0);
    // No atmosphere arithmetic occurs in prepare_wavelength. Both providers
    // size the same worker arrays, and only the factored one restores weights.
    pair.captured.prepare_wavelength(0, 0);
    REQUIRE(std::fegetround() == mode);
    pair.original.prepare_wavelength(0, 0);
    REQUIRE(pair.captured.workspace_bytes() == pair.original.workspace_bytes());
    REQUIRE(std::fesetround(FE_TONEAREST) == 0);
    for (int wavelength = 0; wavelength < wavelengths; ++wavelength) {
        require_products_equal(products(pair.captured, pair.captured_geometry,
                                        atmosphere, wavelength),
                               products(pair.original, pair.original_geometry,
                                        atmosphere, wavelength));
    }
    // The permanent four-double fallback also preserves caller-mode products.
    REQUIRE(std::fesetround(mode) == 0);
    for (int wavelength = 0; wavelength < wavelengths; ++wavelength) {
        require_products_equal(products(pair.captured, pair.captured_geometry,
                                        atmosphere, wavelength),
                               products(pair.original, pair.original_geometry,
                                        atmosphere, wavelength));
        REQUIRE(std::fegetround() == mode);
    }
}

TEST_CASE("Multiworker first-order configurations do not request factored "
          "endpoint storage",
          "[successive_orders][first_order][endpoint_factors]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    FirstOrderProvider<1> provider(geo, tracer);
    auto config = configuration(2);
    provider.initialize_config(config);
    REQUIRE_FALSE(provider.requests_endpoint_factors());
    config.set_threading_model(sasktran2::Config::ThreadingModel::source);
    provider.initialize_config(config);
    REQUIRE_FALSE(provider.requests_endpoint_factors());
}

#endif
