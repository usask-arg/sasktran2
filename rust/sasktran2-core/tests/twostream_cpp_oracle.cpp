#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <sasktran2.h>
#include <sktran_disco/twostream/twostream.h>
#include <string>
#include <string_view>

template <sasktran2::twostream::SourceType Mode>
double run(int iterations = 1) {
    constexpr int n = 3;
    sasktran2::twostream::Input<Mode> input;
    input.nlyr = n;
    input.csz = 0.5;
    input.mu = 0.5;
    input.albedo = 0.2;
    input.thermal_surf = 1.5;
    input.od.resize(n);
    input.ssa.resize(n);
    input.b1.resize(n);
    input.od << 0.01, 0.02, 0.03;
    input.ssa << 0.8, 0.75, 0.7;
    input.b1 << 0.4, 0.3, 0.2;
    if constexpr (sasktran2::twostream::has_solar<Mode>()) {
        input.transmission.resize(n + 1);
        input.average_secant.resize(n);
        input.expsec.resize(n);
        input.average_secant.setConstant(2.0);
        input.transmission << 1.0, std::exp(-0.02), std::exp(-0.06),
            std::exp(-0.12);
        input.expsec = (-input.average_secant.array() * input.od.array()).exp();
    }
    if constexpr (sasktran2::twostream::has_thermal<Mode>()) {
        input.b0_thermal.resize(n);
        input.b1_thermal.resize(n);
        input.exp_thermal.resize(n);
        input.b0_thermal << 2.0, 1.8, 1.6;
        input.b1_thermal << 0.1, 0.15, 0.2;
        input.exp_thermal =
            (-input.b1_thermal.array() * input.od.array()).exp();
    }

    sasktran2::twostream::Solution<Mode> solution;
    sasktran2::twostream::Sources<Mode> sources;
    solution.init(n);
    sources.init(n);
    for (int iteration = 0; iteration < iterations; ++iteration) {
        sasktran2::twostream::solve(input, solution);
        sasktran2::twostream::post_process(input, 0.6, 0.4, solution, sources);
    }

    Eigen::VectorXd weights(n);
    weights(0) = 1.0;
    weights(1) = std::exp(-input.od(0) / 0.6);
    weights(2) = std::exp(-(input.od(0) + input.od(1)) / 0.6);
    double result = sources.source.value.dot(weights);
    double ground = (solution.particular[0].G_plus_bottom.value(n - 1) +
                     solution.bvp_coeffs[0].rhs(2 * n - 2, 0) *
                         solution.homog[0].X_plus.value(n - 1) *
                         solution.homog[0].omega.value(n - 1) +
                     solution.bvp_coeffs[0].rhs(2 * n - 1, 0) *
                         solution.homog[0].X_minus.value(n - 1)) *
                    2 * input.mu * input.albedo * weights(n - 1) *
                    std::exp(-input.od(n - 1) / 0.6);
    if constexpr (sasktran2::twostream::has_thermal<Mode>()) {
        ground += input.thermal_surf * weights(n - 1) *
                  std::exp(-input.od(n - 1) / 0.6);
    }
    return result + ground;
}

double benchmark_vjp(int nlyr, int iterations) {
    constexpr auto mode = sasktran2::twostream::SourceType::ONLY_SOLAR;
    constexpr double solar_cosine = 0.6;
    constexpr double viewing_cosine = 0.6;
    constexpr double azimuth = 0.3;

    sasktran2::Coordinates coords(solar_cosine, 0, 6371000,
                                  sasktran2::geometrytype::planeparallel);
    Eigen::VectorXd grid_values(nlyr + 1);
    for (int i = 0; i <= nlyr; ++i) {
        grid_values(i) = i;
    }
    sasktran2::grids::AltitudeGrid grid(std::move(grid_values),
                                        sasktran2::grids::gridspacing::constant,
                                        sasktran2::grids::outofbounds::extend,
                                        sasktran2::grids::interpolation::lower);
    sasktran2::Geometry1D geometry(std::move(coords), std::move(grid));

    sasktran2::atmosphere::AtmosphereGridStorageFull<1> storage(
        1, geometry.size(), 16);
    sasktran2::atmosphere::Surface<1> surface(1);
    surface.brdf_args().setConstant(0.3 / EIGEN_PI);
    sasktran2::atmosphere::Atmosphere<1> atmosphere(std::move(storage),
                                                    std::move(surface), true);
    atmosphere.storage().total_extinction.setConstant(0.01);
    atmosphere.storage().ssa.setConstant(0.8);
    atmosphere.storage().leg_coeff.chip(0, 0).setConstant(1.0);
    atmosphere.storage().leg_coeff.chip(1, 0).setConstant(0.5);
    atmosphere.storage().solar_irradiance.setConstant(1.0);

    sasktran2::Config config;
    sasktran_disco::PersistentConfiguration<1> persistent;
    std::vector<sasktran2::raytracing::TracedRay> traced_rays;
    sasktran_disco::SKTRAN_DO_UserSpec specs;
    persistent.configure(specs, config, solar_cosine, nlyr, traced_rays);
    sasktran_disco::GeometryLayerArray<1> layer_geometry(persistent, geometry);

    sasktran2::twostream::Input<mode> input;
    input.atmosphere = &atmosphere;
    input.geometry_layers = &layer_geometry;
    input.csz = solar_cosine;
    input.init(nlyr);
    input.calculate(0);

    sasktran2::twostream::Solution<mode> solution;
    sasktran2::twostream::Sources<mode> sources;
    solution.init(nlyr);
    sources.init(nlyr);

    Eigen::MatrixXd attenuation(nlyr, nlyr);
    attenuation.setZero();
    for (int row = 1; row < nlyr; ++row) {
        attenuation(row, Eigen::seq(0, row - 1))
            .setConstant(1.0 / viewing_cosine);
    }
    Eigen::RowVectorXd internal_gradient(3 * (nlyr + 1));
    std::array<Eigen::MatrixXd, 2> bvp_storage;
    for (auto& matrix : bvp_storage) {
        matrix.resize(2 * nlyr, 1);
    }
    sasktran2::Dual<double, sasktran2::dualstorage::dense, 1> output(
        1, atmosphere.num_deriv(), true);

    double checksum = 0.0;
    for (int iteration = 0; iteration < iterations; ++iteration) {
        internal_gradient.setZero();
        output.deriv.setZero();
        sasktran2::twostream::solve(input, solution);
        sources.final_weight_factors.noalias() = (-attenuation * input.od);
        sources.final_weight_factors =
            sources.final_weight_factors.array().exp();
        sasktran2::twostream::post_process(input, viewing_cosine, azimuth,
                                           solution, sources);

        const double ground_weight =
            2 * input.mu * input.albedo *
            sources.final_weight_factors(Eigen::placeholders::last) *
            sources.beamtrans.value(Eigen::placeholders::last);
        const double ground_source =
            (solution.particular[0].G_plus_bottom.value(
                 Eigen::placeholders::last) +
             solution.bvp_coeffs[0].rhs(Eigen::placeholders::last - 1, 0) *
                 solution.homog[0].X_plus.value(Eigen::placeholders::last) *
                 solution.homog[0].omega.value(Eigen::placeholders::last) +
             solution.bvp_coeffs[0].rhs(Eigen::placeholders::last, 0) *
                 solution.homog[0].X_minus.value(Eigen::placeholders::last)) *
            ground_weight;

        sasktran2::twostream::backprop::GradientMap<mode> gradient(
            atmosphere, internal_gradient.data());
        gradient.d_albedo = ground_source / input.albedo;
        sasktran2::twostream::backprop::full(
            input, solution, sources, sources.final_weight_factors, bvp_storage,
            gradient, ground_weight);
        for (int layer = 0; layer < nlyr; ++layer) {
            sasktran2::twostream::backprop::od(
                input,
                -sources.final_weight_factors(layer) *
                    sources.source.value(layer) * attenuation.row(layer),
                gradient);
        }
        sasktran2::twostream::backprop::od(
            input,
            -ground_source * Eigen::RowVectorXd::Ones(nlyr) / viewing_cosine,
            gradient);
        sasktran2::twostream::backprop::map_to_atmosphere(input, gradient,
                                                          output);
        checksum += output.deriv(0);
    }
    return checksum;
}

struct ThreadBenchmark {
    double checksum;
    double nanoseconds_per_batch;
};

ThreadBenchmark benchmark_vjp_threads(int nlyr, int wavelengths, int batches,
                                      int threads) {
    double checksum = 0.0;
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    const auto start = std::chrono::steady_clock::now();
#ifdef _OPENMP
#pragma omp parallel num_threads(threads) reduction(+ : checksum)
    {
        const int thread = omp_get_thread_num();
        const int team_size = omp_get_num_threads();
        const int wavelengths_for_thread =
            wavelengths / team_size + (thread < wavelengths % team_size);
        checksum += benchmark_vjp(nlyr, wavelengths_for_thread * batches);
    }
#else
    checksum = benchmark_vjp(nlyr, wavelengths * batches);
#endif
    const auto elapsed = std::chrono::duration<double, std::nano>(
        std::chrono::steady_clock::now() - start);
    return {checksum, elapsed.count() / batches};
}

int main(int argc, char** argv) {
    if (argc > 1) {
        if (std::string_view(argv[1]) == "vjp_mt") {
            const int nlyr = argc > 2 ? std::stoi(argv[2]) : 20;
            const int wavelengths = argc > 3 ? std::stoi(argv[3]) : 8192;
            const int batches = argc > 4 ? std::stoi(argv[4]) : 20;
            const int threads = argc > 5 ? std::stoi(argv[5]) : 1;
            const auto result =
                benchmark_vjp_threads(nlyr, wavelengths, batches, threads);
            std::cout << result.checksum << "\n"
                      << result.nanoseconds_per_batch << " ns/batch\n";
            return 0;
        }
        const bool derivatives = std::string_view(argv[1]) == "vjp";
        const int nlyr = argc > 2 ? std::stoi(argv[2]) : 3;
        const int iterations = argc > 3 ? std::stoi(argv[3]) : 100000;
        const auto start = std::chrono::steady_clock::now();
        volatile double result =
            derivatives
                ? benchmark_vjp(nlyr, iterations)
                : run<sasktran2::twostream::SourceType::ONLY_SOLAR>(iterations);
        const auto elapsed = std::chrono::duration<double, std::nano>(
            std::chrono::steady_clock::now() - start);
        std::cout << result << "\n"
                  << elapsed.count() / iterations << " ns/solve\n";
        return 0;
    }
    std::cout << std::setprecision(17)
              << run<sasktran2::twostream::SourceType::ONLY_SOLAR>() << "\n"
              << run<sasktran2::twostream::SourceType::ONLY_THERMAL>() << "\n";
}
