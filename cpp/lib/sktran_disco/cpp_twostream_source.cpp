#include "sktran_disco/twostream/cpp_source.h"

#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_geometrylayerarray.h"
#include "sktran_disco/sktran_do_pconfig.h"
#include "sktran_disco/sktran_do_specs.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace {
    constexpr int LANES = 8;
    constexpr double FOUR_PI = 4.0 * EIGEN_PI;
    constexpr double THERMAL_MIN_OPTICAL_DEPTH = 1.0e-10;
    constexpr double THERMAL_MIN_EMISSION = 1.0e-30;

    using Wide = Eigen::Array<double, LANES, 1>;

    struct StableValue {
        double value;
        double da;
        double db;
        double dt;
    };

    double exp_moment(int order, double rate, double thickness) {
        if (thickness == 0.0) {
            return 0.0;
        }
        const double scaled_rate = rate * thickness;
        double unit_moment;
        if (std::abs(scaled_rate) < 0.5) {
            double factorial_term = 1.0;
            unit_moment = 0.0;
            for (int term = 0; term < 40; ++term) {
                const double contribution =
                    factorial_term / static_cast<double>(order + term + 1);
                unit_moment += contribution;
                if (std::abs(contribution) <=
                    std::numeric_limits<double>::epsilon() *
                        std::max(std::abs(unit_moment), 1.0)) {
                    break;
                }
                factorial_term *= -scaled_rate / static_cast<double>(term + 1);
            }
        } else {
            const double exponential = std::exp(-scaled_rate);
            unit_moment = -std::expm1(-scaled_rate) / scaled_rate;
            for (int current = 1; current <= order; ++current) {
                unit_moment =
                    (static_cast<double>(current) * unit_moment - exponential) /
                    scaled_rate;
            }
        }
        return std::pow(thickness, order + 1) * unit_moment;
    }

    StableValue exp_difference(double a, double b, double thickness) {
        const double delta = b - a;
        const double scaled_delta = delta * thickness;
        if (std::abs(scaled_delta) > 1.0e-5) {
            const double exp_a = std::exp(-a * thickness);
            const double exp_b = std::exp(-b * thickness);
            const double numerator = exp_a - exp_b;
            const double inv_delta = 1.0 / delta;
            const double inv_delta_squared = inv_delta * inv_delta;
            return {numerator * inv_delta,
                    (-thickness * exp_a * delta + numerator) *
                        inv_delta_squared,
                    (thickness * exp_b * delta - numerator) * inv_delta_squared,
                    (-a * exp_a + b * exp_b) * inv_delta};
        }

        const double midpoint = 0.5 * (a + b);
        const double half_delta = 0.5 * delta;
        const double u = half_delta * thickness;
        const double u2 = u * u;
        const double sinhc =
            1.0 + u2 * (1.0 / 6.0 + u2 * (1.0 / 120.0 + u2 / 5040.0));
        const double sinhc_prime =
            u * (1.0 / 3.0 + u2 * (1.0 / 30.0 + u2 / 840.0));
        const double exp_midpoint = std::exp(-midpoint * thickness);
        const double value = thickness * exp_midpoint * sinhc;
        const double d_midpoint = -thickness * value;
        const double d_half_delta =
            thickness * thickness * exp_midpoint * sinhc_prime;
        return {value, 0.5 * (d_midpoint - d_half_delta),
                0.5 * (d_midpoint + d_half_delta),
                exp_midpoint * ((1.0 - midpoint * thickness) * sinhc +
                                thickness * half_delta * sinhc_prime)};
    }

    StableValue integrated_exp_difference(double a, double b,
                                          double thickness) {
        const double delta = b - a;
        if (std::abs(delta * thickness) > 1.0e-4) {
            const StableValue ga = exp_difference(0.0, a, thickness);
            const StableValue gb = exp_difference(0.0, b, thickness);
            const double numerator = ga.value - gb.value;
            const double inv_delta = 1.0 / delta;
            const double inv_delta_squared = inv_delta * inv_delta;
            return {numerator * inv_delta,
                    (ga.db * delta + numerator) * inv_delta_squared,
                    (-gb.db * delta - numerator) * inv_delta_squared,
                    exp_difference(a, b, thickness).value};
        }

        const double midpoint = 0.5 * (a + b);
        const double half_delta = 0.5 * delta;
        const double half_delta_squared = half_delta * half_delta;
        const double moment1 = exp_moment(1, midpoint, thickness);
        const double moment2 = exp_moment(2, midpoint, thickness);
        const double moment3 = exp_moment(3, midpoint, thickness);
        const double moment4 = exp_moment(4, midpoint, thickness);
        const double value = moment1 + half_delta_squared * moment3 / 6.0;
        const double d_midpoint = -moment2 - half_delta_squared * moment4 / 6.0;
        const double d_half_delta = half_delta * moment3 / 3.0;
        return {value, 0.5 * (d_midpoint - d_half_delta),
                0.5 * (d_midpoint + d_half_delta),
                exp_difference(a, b, thickness).value};
    }

    struct Tape;

    struct Var {
        std::size_t index = 0;
        Tape* tape = nullptr;

        const Wide& value() const;
        Var exp() const;
        Var sqrt() const;
    };

    enum class Op : std::uint8_t {
        Leaf,
        Add,
        Subtract,
        Multiply,
        Divide,
        Negate,
        Scale,
        Exp,
        Sqrt,
        ConstantDivide,
        ExpDifference,
        IntegratedExpDifference,
        PositiveRatio,
        ClampSsa,
        ThermalSlope,
    };

    struct Node {
        std::array<std::size_t, 3> input{};
        double parameter = 0.0;
        Op op = Op::Leaf;
    };

    struct Tape {
        std::vector<Wide> values;
        std::vector<Wide> adjoints;
        std::vector<Node> nodes;
        bool recording = true;

        void reset(std::size_t capacity, bool record_derivatives) {
            values.clear();
            adjoints.clear();
            nodes.clear();
            recording = record_derivatives;
            values.reserve(capacity);
            if (recording) {
                adjoints.reserve(capacity);
                nodes.reserve(capacity);
            }
        }

        Var push(const Wide& value, Node node = {}) {
            const std::size_t index = values.size();
            values.push_back(value);
            if (recording) {
                adjoints.emplace_back(Wide::Zero());
                nodes.push_back(std::move(node));
            }
            return {index, this};
        }

        Var leaf(const Wide& value) { return push(value); }
        Var constant(double value) { return leaf(Wide::Constant(value)); }

        Var unary(const Var& input, const Wide& value, Op op,
                  double parameter = 0.0) {
            if (!recording) {
                return push(value);
            }
            Node node;
            node.input[0] = input.index;
            node.parameter = parameter;
            node.op = op;
            return push(value, std::move(node));
        }

        Var binary(const Var& left, const Var& right, const Wide& value,
                   Op op) {
            if (!recording) {
                return push(value);
            }
            Node node;
            node.input[0] = left.index;
            node.input[1] = right.index;
            node.op = op;
            return push(value, std::move(node));
        }

        Var ternary(const Var& a, const Var& b, const Var& c, const Wide& value,
                    Op op) {
            if (!recording) {
                return push(value);
            }
            Node node;
            node.input = {a.index, b.index, c.index};
            node.op = op;
            return push(value, std::move(node));
        }

        void reverse(std::size_t output) {
            std::fill(adjoints.begin(), adjoints.end(), Wide::Zero());
            adjoints[output].setOnes();
            for (std::size_t index = nodes.size(); index-- > 0;) {
                const Node& node = nodes[index];
                const Wide seed = adjoints[index];
                const std::size_t a = node.input[0];
                const std::size_t b = node.input[1];
                const std::size_t c = node.input[2];
                switch (node.op) {
                case Op::Leaf:
                    break;
                case Op::Add:
                    adjoints[a] += seed;
                    adjoints[b] += seed;
                    break;
                case Op::Subtract:
                    adjoints[a] += seed;
                    adjoints[b] -= seed;
                    break;
                case Op::Multiply:
                    adjoints[a] += seed * values[b];
                    adjoints[b] += seed * values[a];
                    break;
                case Op::Divide: {
                    const Wide inverse = 1.0 / values[b];
                    adjoints[a] += seed * inverse;
                    adjoints[b] -= seed * values[a] * inverse.square();
                    break;
                }
                case Op::Negate:
                    adjoints[a] -= seed;
                    break;
                case Op::Scale:
                    adjoints[a] += seed * node.parameter;
                    break;
                case Op::Exp:
                    adjoints[a] += seed * values[index];
                    break;
                case Op::Sqrt:
                    adjoints[a] += seed * 0.5 / values[index];
                    break;
                case Op::ConstantDivide:
                    adjoints[a] -= seed * node.parameter / values[a].square();
                    break;
                case Op::ExpDifference:
                case Op::IntegratedExpDifference: {
                    Wide da;
                    Wide db;
                    Wide dt;
                    for (int lane = 0; lane < LANES; ++lane) {
                        const StableValue result =
                            node.op == Op::ExpDifference
                                ? exp_difference(values[a][lane],
                                                 values[b][lane],
                                                 values[c][lane])
                                : integrated_exp_difference(values[a][lane],
                                                            values[b][lane],
                                                            values[c][lane]);
                        da[lane] = result.da;
                        db[lane] = result.db;
                        dt[lane] = result.dt;
                    }
                    adjoints[a] += seed * da;
                    adjoints[b] += seed * db;
                    adjoints[c] += seed * dt;
                    break;
                }
                case Op::PositiveRatio:
                    for (int lane = 0; lane < LANES; ++lane) {
                        if (values[b][lane] > 0.0) {
                            adjoints[a][lane] += seed[lane] / values[b][lane];
                            adjoints[b][lane] -=
                                seed[lane] * values[a][lane] /
                                (values[b][lane] * values[b][lane]);
                        }
                    }
                    break;
                case Op::ClampSsa:
                    for (int lane = 0; lane < LANES; ++lane) {
                        if (values[a][lane] < 1.0 - 1.0e-9) {
                            adjoints[a][lane] += seed[lane];
                        }
                    }
                    break;
                case Op::ThermalSlope:
                    for (int lane = 0; lane < LANES; ++lane) {
                        const double top = values[a][lane];
                        const double bottom = values[b][lane];
                        const double od = values[c][lane];
                        if (od > THERMAL_MIN_OPTICAL_DEPTH &&
                            top > THERMAL_MIN_EMISSION &&
                            bottom > THERMAL_MIN_EMISSION) {
                            adjoints[a][lane] += seed[lane] / (top * od);
                            adjoints[b][lane] -= seed[lane] / (bottom * od);
                            adjoints[c][lane] -=
                                seed[lane] * values[index][lane] / od;
                        }
                    }
                    break;
                }
            }
        }
    };

    const Wide& Var::value() const { return tape->values[index]; }

    Var Var::exp() const {
        const Wide result = value().exp();
        return tape->unary(*this, result, Op::Exp);
    }

    Var Var::sqrt() const {
        const Wide result = value().sqrt();
        if (!tape->recording) {
            return tape->push(result);
        }
        return tape->unary(*this, result, Op::Sqrt);
    }

    Var operator+(const Var& a, const Var& b) {
        return a.tape->binary(a, b, a.value() + b.value(), Op::Add);
    }
    Var operator-(const Var& a, const Var& b) {
        return a.tape->binary(a, b, a.value() - b.value(), Op::Subtract);
    }
    Var operator*(const Var& a, const Var& b) {
        return a.tape->binary(a, b, a.value() * b.value(), Op::Multiply);
    }
    Var operator/(const Var& a, const Var& b) {
        const Wide inverse = 1.0 / b.value();
        if (!a.tape->recording) {
            return a.tape->push(a.value() * inverse);
        }
        return a.tape->binary(a, b, a.value() * inverse, Op::Divide);
    }
    Var operator-(const Var& value) {
        return value.tape->unary(value, -value.value(), Op::Negate);
    }

    Var affine(const Var& value, double scale, double offset) {
        return value.tape->unary(value, value.value() * scale + offset,
                                 Op::Scale, scale);
    }
    Var operator+(const Var& a, double b) { return affine(a, 1.0, b); }
    Var operator+(double a, const Var& b) { return b + a; }
    Var operator-(const Var& a, double b) { return affine(a, 1.0, -b); }
    Var operator-(double a, const Var& b) { return affine(b, -1.0, a); }
    Var operator*(const Var& a, double b) { return affine(a, b, 0.0); }
    Var operator*(double a, const Var& b) { return b * a; }
    Var operator/(const Var& a, double b) { return a * (1.0 / b); }
    Var operator/(double a, const Var& b) {
        const Wide inverse = 1.0 / b.value();
        if (!b.tape->recording) {
            return b.tape->push(a * inverse);
        }
        return b.tape->unary(b, a * inverse, Op::ConstantDivide, a);
    }

    Var stable_exp_difference(const Var& a, const Var& b, const Var& t) {
        Wide value;
        for (int lane = 0; lane < LANES; ++lane) {
            const StableValue result = exp_difference(
                a.value()[lane], b.value()[lane], t.value()[lane]);
            value[lane] = result.value;
        }
        return a.tape->ternary(a, b, t, value, Op::ExpDifference);
    }

    Var stable_integrated_exp_difference(const Var& a, const Var& b,
                                         const Var& t) {
        Wide value;
        for (int lane = 0; lane < LANES; ++lane) {
            const StableValue result = integrated_exp_difference(
                a.value()[lane], b.value()[lane], t.value()[lane]);
            value[lane] = result.value;
        }
        return a.tape->ternary(a, b, t, value, Op::IntegratedExpDifference);
    }

    Var positive_ratio(const Var& numerator, const Var& denominator) {
        Wide value = Wide::Zero();
        for (int lane = 0; lane < LANES; ++lane) {
            if (denominator.value()[lane] > 0.0) {
                value[lane] =
                    numerator.value()[lane] / denominator.value()[lane];
            }
        }
        return numerator.tape->binary(numerator, denominator, value,
                                      Op::PositiveRatio);
    }

    Var clamp_ssa(const Var& value) {
        constexpr double upper = 1.0 - 1.0e-9;
        Wide result;
        for (int lane = 0; lane < LANES; ++lane) {
            result[lane] = std::min(value.value()[lane], upper);
        }
        return value.tape->unary(value, result, Op::ClampSsa);
    }

    Var thermal_profile_slope(const Var& top, const Var& bottom,
                              const Var& optical_depth) {
        Wide value = Wide::Zero();
        for (int lane = 0; lane < LANES; ++lane) {
            const double t = top.value()[lane];
            const double b = bottom.value()[lane];
            const double od = optical_depth.value()[lane];
            if (od > THERMAL_MIN_OPTICAL_DEPTH && t > THERMAL_MIN_EMISSION &&
                b > THERMAL_MIN_EMISSION) {
                value[lane] = std::log1p((t - b) / b) / od;
            }
        }
        return top.tape->ternary(top, bottom, optical_depth, value,
                                 Op::ThermalSlope);
    }

    void resize_vars(std::vector<Var>& values, std::size_t size,
                     const Var& zero) {
        values.resize(size, zero);
        std::fill(values.begin(), values.end(), zero);
    }

    struct ColumnGeometry {
        std::vector<double> layer_thickness;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
            chapman;
        double solar_cosine = 0.0;
        double quadrature_cosine = 0.5;
    };

    struct ViewGeometry {
        double cosine = 0.0;
        double relative_azimuth = 0.0;
    };

    struct SphericalGeometry {
        std::vector<double> sza_grid;
        std::vector<ColumnGeometry> columns;
        std::vector<std::size_t> ray_offsets;
        std::vector<std::uint8_t> ground_hit;
        std::vector<double> ground_cos_sza;
        std::vector<std::size_t> segment_layers;
        std::vector<double> segment_fractions;
        std::vector<double> segment_cosines;
        std::vector<double> segment_relative_azimuths;
        std::vector<double> segment_cos_sza;
        std::vector<std::size_t> od_offsets;
        std::vector<std::size_t> od_indices;
        std::vector<double> od_weights;
    };

    struct Homogeneous {
        std::vector<Var> k;
        std::vector<Var> xp;
        std::vector<Var> xm;
        std::vector<Var> omega;
    };

    struct Particular {
        std::vector<Var> ap;
        std::vector<Var> am;
        std::vector<Var> at;
        std::vector<Var> gpt;
        std::vector<Var> gpb;
        std::vector<Var> gmt;
        std::vector<Var> gmb;
    };

    struct Bvp {
        std::vector<Var> e, c, d, a, b, rhs;
        std::vector<Var> alpha, beta, gamma, inverse_mu, z;
    };

    struct ColumnSolution {
        std::vector<Var> od;
        std::vector<Var> ssa;
        std::vector<Var> b1;
        std::vector<Var> transmission;
        std::vector<Var> secant;
        std::vector<Var> thermal_b0;
        std::vector<Var> thermal_b1;
        std::array<Homogeneous, 2> homogeneous;
        std::array<Particular, 2> particular;
        std::array<Bvp, 2> bvp;
    };

    struct PacketWorkspace {
        Tape tape;
        std::vector<Var> extinction;
        std::vector<Var> level_ssa;
        std::vector<Var> level_b1;
        std::vector<Var> emission;
        std::vector<ColumnSolution> columns;
        std::vector<Var> outputs;
    };

    Wide positive_ratio(const Wide& numerator, const Wide& denominator) {
        Wide result = Wide::Zero();
        for (int lane = 0; lane < LANES; ++lane) {
            if (denominator[lane] > 0.0) {
                result[lane] = numerator[lane] / denominator[lane];
            }
        }
        return result;
    }

    Wide positive_reciprocal(const Wide& value) {
        Wide result = Wide::Zero();
        for (int lane = 0; lane < LANES; ++lane) {
            if (value[lane] > 0.0) {
                result[lane] = 1.0 / value[lane];
            }
        }
        return result;
    }

    Wide stable_exp_difference(const Wide& a, const Wide& b,
                               const Wide& thickness) {
        Wide result;
        for (int lane = 0; lane < LANES; ++lane) {
            result[lane] =
                exp_difference(a[lane], b[lane], thickness[lane]).value;
        }
        return result;
    }

    Wide stable_integrated_exp_difference(const Wide& a, const Wide& b,
                                          const Wide& thickness) {
        Wide result;
        for (int lane = 0; lane < LANES; ++lane) {
            result[lane] =
                integrated_exp_difference(a[lane], b[lane], thickness[lane])
                    .value;
        }
        return result;
    }

    std::tuple<Wide, Wide, Wide> exp_difference_adjoint(const Wide& a,
                                                        const Wide& b,
                                                        const Wide& thickness,
                                                        const Wide& seed) {
        Wide da;
        Wide db;
        Wide dt;
        for (int lane = 0; lane < LANES; ++lane) {
            const StableValue result =
                exp_difference(a[lane], b[lane], thickness[lane]);
            da[lane] = seed[lane] * result.da;
            db[lane] = seed[lane] * result.db;
            dt[lane] = seed[lane] * result.dt;
        }
        return {da, db, dt};
    }

    bool exp_difference_ratio_nonresonant(const Wide& a, const Wide& b,
                                          const Wide& thickness) {
        for (int lane = 0; lane < LANES; ++lane) {
            if (std::abs((b[lane] - a[lane]) * thickness[lane]) <= 1.0e-5) {
                return false;
            }
        }
        return true;
    }

    Wide exp_difference_ratio(const Wide& exp_a, const Wide& exp_b,
                              const Wide& a, const Wide& b,
                              const Wide& thickness) {
        if (exp_difference_ratio_nonresonant(a, b, thickness)) {
            return (exp_a - exp_b) / (b - a);
        }
        return stable_exp_difference(a, b, thickness);
    }

    std::tuple<Wide, Wide, Wide>
    exp_difference_ratio_adjoint(const Wide& exp_a, const Wide& exp_b,
                                 const Wide& a, const Wide& b,
                                 const Wide& thickness, const Wide& seed) {
        if (!exp_difference_ratio_nonresonant(a, b, thickness)) {
            return exp_difference_adjoint(a, b, thickness, seed);
        }
        const Wide delta = b - a;
        const Wide numerator = exp_a - exp_b;
        const Wide inverse_delta = 1.0 / delta;
        const Wide inverse_delta_squared = inverse_delta.square();
        return {seed * (-thickness * exp_a * delta + numerator) *
                    inverse_delta_squared,
                seed * (thickness * exp_b * delta - numerator) *
                    inverse_delta_squared,
                seed * (-a * exp_a + b * exp_b) * inverse_delta};
    }

    std::tuple<Wide, Wide, Wide>
    integrated_exp_difference_adjoint(const Wide& a, const Wide& b,
                                      const Wide& thickness, const Wide& seed) {
        Wide da;
        Wide db;
        Wide dt;
        for (int lane = 0; lane < LANES; ++lane) {
            const StableValue result =
                integrated_exp_difference(a[lane], b[lane], thickness[lane]);
            da[lane] = seed[lane] * result.da;
            db[lane] = seed[lane] * result.db;
            dt[lane] = seed[lane] * result.dt;
        }
        return {da, db, dt};
    }

    Wide thermal_profile_slope(const Wide& top, const Wide& bottom,
                               const Wide& optical_depth) {
        Wide result = Wide::Zero();
        for (int lane = 0; lane < LANES; ++lane) {
            if (optical_depth[lane] > THERMAL_MIN_OPTICAL_DEPTH &&
                top[lane] > THERMAL_MIN_EMISSION &&
                bottom[lane] > THERMAL_MIN_EMISSION) {
                result[lane] =
                    std::log1p((top[lane] - bottom[lane]) / bottom[lane]) /
                    optical_depth[lane];
            }
        }
        return result;
    }

    struct ExplicitHomogeneous {
        std::vector<Wide> d, s, k, xp, xm, omega, norm;
    };

    struct ExplicitHomogeneousAdjoint {
        std::vector<Wide> k, xp, xm, omega;
    };

    struct ExplicitParticular {
        std::vector<Wide> qp, qm, ap, am, at, exponential, cp, cm;
        std::vector<Wide> gpt, gpb, gmt, gmb;
    };

    struct ExplicitParticularAdjoint {
        std::vector<Wide> ap, am, at, gpt, gpb, gmt, gmb;
    };

    struct ExplicitBvp {
        std::vector<Wide> e, c, d, a, b, rhs;
        std::vector<Wide> inverse_mu, alpha, beta, gamma, z;
    };

    struct ExplicitWorkspace {
        std::vector<Wide> level_extinction, level_ssa, level_b1, level_emission;
        std::vector<Wide> od, ssa, b1, transmission, secant;
        std::vector<Wide> thermal_b0, thermal_b1;
        std::vector<Wide> d_od, d_ssa, d_b1, d_transmission, d_secant;
        std::vector<Wide> d_thermal_b0, d_thermal_b1;
        std::vector<Wide> d_level_extinction, d_level_ssa, d_level_b1,
            d_level_emission;
        std::vector<Wide> beam, source, attenuation;
        std::array<ExplicitHomogeneous, 2> homogeneous;
        std::array<ExplicitHomogeneousAdjoint, 2> d_homogeneous;
        std::array<ExplicitParticular, 2> particular;
        std::array<ExplicitParticularAdjoint, 2> d_particular;
        std::array<ExplicitBvp, 2> bvp;
        std::array<std::vector<Wide>, 2> d_solution;
    };

    struct SphericalExplicitWorkspace {
        std::vector<ExplicitWorkspace> columns;
        std::vector<Wide> path_radiance;
        std::vector<Wide> path_transmission;
        std::vector<Wide> path_source;
        std::vector<Wide> d_level_extinction;
        std::vector<Wide> d_level_ssa;
        std::vector<Wide> d_level_b1;
        std::vector<Wide> d_level_emission;
    };

    void resize_wide(std::vector<Wide>& values, int size, bool zero);

    void resize_explicit_workspace(ExplicitWorkspace& workspace, int n,
                                   bool solar, bool clear_adjoints);

    template <bool Solar>
    void prepare_explicit_column(const ColumnGeometry& geometry,
                                 const Wide& irradiance,
                                 ExplicitWorkspace& workspace);

    template <bool Solar>
    void forward_explicit_layers(const ColumnGeometry& geometry,
                                 const Wide& albedo,
                                 const Wide& thermal_surface,
                                 ExplicitWorkspace& workspace);

    template <bool Solar>
    bool explicit_packet_nonresonant(const ColumnGeometry& geometry,
                                     const std::vector<ViewGeometry>& views,
                                     const ExplicitWorkspace& workspace);

    template <bool Solar>
    Wide explicit_plane_view(const ColumnGeometry& geometry,
                             const ViewGeometry& view, const Wide& albedo,
                             const Wide& thermal_surface, bool reverse,
                             Wide& d_albedo, Wide& d_thermal_surface,
                             ExplicitWorkspace& workspace);

    template <bool Solar>
    std::pair<Wide, Wide> reverse_explicit_bvp(const ColumnGeometry& geometry,
                                               const Wide& albedo,
                                               ExplicitWorkspace& workspace);

    template <bool Solar>
    void reverse_explicit_layers(const ColumnGeometry& geometry,
                                 ExplicitWorkspace& workspace);

    template <bool Solar>
    void map_explicit_to_levels(const ColumnGeometry& geometry,
                                ExplicitWorkspace& workspace);

    template <bool Solar>
    bool explicit_column_nonresonant(const ColumnGeometry& geometry,
                                     const ExplicitWorkspace& workspace);

    template <bool Solar>
    Wide explicit_local_source(const ColumnGeometry& geometry, int layer,
                               double fraction_from_top,
                               double propagation_cosine,
                               double relative_azimuth,
                               const ExplicitWorkspace& workspace);

    template <bool Solar>
    void reverse_explicit_local_source(const ColumnGeometry& geometry,
                                       int layer, double fraction_from_top,
                                       double propagation_cosine,
                                       double relative_azimuth,
                                       const Wide& seed,
                                       ExplicitWorkspace& workspace);

    template <bool Solar>
    Wide explicit_surface_source(const ColumnGeometry& geometry,
                                 const Wide& albedo,
                                 const Wide& thermal_surface,
                                 const ExplicitWorkspace& workspace);

    template <bool Solar>
    std::pair<Wide, Wide>
    reverse_explicit_surface_source(const ColumnGeometry& geometry,
                                    const Wide& albedo, const Wide& seed,
                                    ExplicitWorkspace& workspace);

    std::pair<double, double>
    normalized_source_weights(const sasktran2::raytracing::TracedLayer& layer) {
        const double sum =
            layer.od_quad_start_fraction + layer.od_quad_end_fraction;
        if (std::isfinite(layer.od_quad_start_fraction) &&
            std::isfinite(layer.od_quad_end_fraction) && sum > 0.0) {
            return {layer.od_quad_start_fraction / sum,
                    layer.od_quad_end_fraction / sum};
        }
        return {0.5, 0.5};
    }

    template <bool Solar>
    void prepare_and_solve_column(const ColumnGeometry& geometry,
                                  const std::vector<Var>& extinction,
                                  const std::vector<Var>& level_ssa,
                                  const std::vector<Var>& level_b1,
                                  const std::vector<Var>& emission,
                                  const Var& irradiance, const Var& albedo,
                                  const Var& thermal_surface,
                                  ColumnSolution& solution, const Var& zero);

    template <bool Solar>
    Var plane_radiance(const ColumnGeometry& geometry,
                       const ColumnSolution& solution, const ViewGeometry& view,
                       const Var& albedo, const Var& thermal_surface,
                       const Var& zero);

    template <bool Solar>
    Var spherical_radiance(const SphericalGeometry& spherical,
                           const std::vector<ColumnSolution>& columns,
                           const std::vector<Var>& extinction, int view,
                           const Var& albedo, const Var& thermal_surface,
                           const Var& zero);
} // namespace

template <sasktran2::twostream::SourceType SOURCE_TYPE>
struct CppTwoStreamSourceAdapter<SOURCE_TYPE>::Impl {
    struct Worker {
        PacketWorkspace packet;
        ExplicitWorkspace explicit_plane;
        SphericalExplicitWorkspace explicit_spherical;
        std::vector<double> radiance;
        std::vector<double> extinction;
        std::vector<double> ssa;
        std::vector<double> b1;
        std::vector<double> emission;
        std::vector<double> surface_albedo;
        std::vector<double> surface_emission;
    };

    const sasktran2::Geometry1D& geometry;
    const sasktran2::Config* config = nullptr;
    const sasktran2::atmosphere::Atmosphere<1>* atmosphere = nullptr;
    sasktran_disco::SKTRAN_DO_UserSpec spec;
    std::vector<std::unique_ptr<sasktran_disco::PersistentConfiguration<1>>>
        pconfigs;
    std::vector<std::unique_ptr<sasktran_disco::GeometryLayerArray<1>>>
        geometry_layers;
    std::vector<ColumnGeometry> columns;
    std::vector<ViewGeometry> views;
    std::unique_ptr<SphericalGeometry> spherical;
    std::vector<Worker> workers;
    int block_capacity = 1;

    explicit Impl(const sasktran2::Geometry1D& geometry) : geometry(geometry) {
        spec.configure(2, geometry.size() - 1);
    }

    static Wide load_matrix(const Eigen::Ref<const Eigen::MatrixXd>& values,
                            int level,
                            const sasktran2::WavelengthBlock<>& block,
                            int packet_base) {
        Wide result;
        for (int lane = 0; lane < LANES; ++lane) {
            const int local = std::min(packet_base + lane, block.count - 1);
            result[lane] = values(level, block.wavelength(local));
        }
        return result;
    }

    Wide load_extinction(int top_down_level,
                         const sasktran2::WavelengthBlock<>& block,
                         int packet_base) const {
        const int cpp_level = geometry.size() - 1 - top_down_level;
        return load_matrix(atmosphere->storage().total_extinction, cpp_level,
                           block, packet_base);
    }

    Wide load_ssa(int top_down_level, const sasktran2::WavelengthBlock<>& block,
                  int packet_base) const {
        const int cpp_level = geometry.size() - 1 - top_down_level;
        return load_matrix(atmosphere->storage().ssa, cpp_level, block,
                           packet_base);
    }

    Wide load_b1(int top_down_level, const sasktran2::WavelengthBlock<>& block,
                 int packet_base) const {
        const int cpp_level = geometry.size() - 1 - top_down_level;
        Wide result;
        for (int lane = 0; lane < LANES; ++lane) {
            const int local = std::min(packet_base + lane, block.count - 1);
            const int wavelength = block.wavelength(local);
            const double delta_m =
                atmosphere->storage().f(cpp_level, wavelength);
            const double denominator = 1.0 - delta_m;
            if (!std::isfinite(delta_m) || denominator == 0.0) {
                throw std::invalid_argument(
                    "delta-M fraction contains an invalid value");
            }
            result[lane] =
                atmosphere->storage().leg_coeff(1, cpp_level, wavelength) -
                3.0 * delta_m / denominator;
        }
        return result;
    }

    Wide load_emission(int top_down_level,
                       const sasktran2::WavelengthBlock<>& block,
                       int packet_base) const {
        const int cpp_level = geometry.size() - 1 - top_down_level;
        return load_matrix(atmosphere->storage().emission_source, cpp_level,
                           block, packet_base);
    }

    Wide load_spectral(const Eigen::Ref<const Eigen::VectorXd>& values,
                       const sasktran2::WavelengthBlock<>& block,
                       int packet_base) const {
        Wide result;
        for (int lane = 0; lane < LANES; ++lane) {
            const int local = std::min(packet_base + lane, block.count - 1);
            result[lane] = values[block.wavelength(local)];
        }
        return result;
    }

    Wide load_albedo(const sasktran2::WavelengthBlock<>& block,
                     int packet_base) const {
        Wide result;
        for (int lane = 0; lane < LANES; ++lane) {
            const int local = std::min(packet_base + lane, block.count - 1);
            result[lane] = atmosphere->surface().brdf(block.wavelength(local),
                                                      0, 0, 0)(0, 0) *
                           EIGEN_PI;
        }
        return result;
    }

    void resize_worker(Worker& worker) const {
        const std::size_t nviews =
            spherical ? spherical->ground_hit.size() : views.size();
        const std::size_t nlevels = geometry.size();
        const std::size_t surface_size =
            nviews * static_cast<std::size_t>(block_capacity);
        const std::size_t level_size = surface_size * nlevels;
        worker.radiance.resize(surface_size);
        worker.surface_albedo.resize(surface_size);
        worker.surface_emission.resize(surface_size);
        worker.extinction.resize(level_size);
        worker.ssa.resize(level_size);
        worker.b1.resize(level_size);
        worker.emission.resize(level_size);
    }

    template <bool Solar>
    bool calculate_explicit_packet(const sasktran2::WavelengthBlock<>& block,
                                   int packet_base, Worker& worker) {
        const int nlevels = geometry.size();
        const int n = nlevels - 1;
        const int nviews = static_cast<int>(views.size());
        const int valid_lanes = std::min(LANES, block.count - packet_base);
        const bool jacobians = atmosphere->num_deriv() > 0;
        auto& workspace = worker.explicit_plane;
        // Forward state is overwritten below.  Adjoint storage is cleared only
        // when a view actually requests a reverse sweep.
        resize_explicit_workspace(workspace, n, Solar, false);
        for (int level = 0; level < nlevels; ++level) {
            workspace.level_extinction[level] =
                load_extinction(level, block, packet_base);
            workspace.level_ssa[level] = load_ssa(level, block, packet_base);
            workspace.level_b1[level] = load_b1(level, block, packet_base);
            if constexpr (!Solar) {
                workspace.level_emission[level] =
                    load_emission(level, block, packet_base);
            }
        }
        const Wide albedo = load_albedo(block, packet_base);
        Wide irradiance = Wide::Zero();
        Wide thermal_surface = Wide::Zero();
        if constexpr (Solar) {
            irradiance = load_spectral(atmosphere->storage().solar_irradiance,
                                       block, packet_base);
        } else {
            thermal_surface = load_spectral(atmosphere->surface().emission(),
                                            block, packet_base);
        }
        prepare_explicit_column<Solar>(columns[0], irradiance, workspace);
        forward_explicit_layers<Solar>(columns[0], albedo, thermal_surface,
                                       workspace);

        // The closed-form reverse kernel uses direct quotient formulas.  The
        // stable operation tape remains the fallback for removable
        // singularities, while the forward-only kernel is stable everywhere.
        if (jacobians &&
            !explicit_packet_nonresonant<Solar>(columns[0], views, workspace)) {
            return false;
        }

        for (int view = 0; view < nviews; ++view) {
            if (jacobians) {
                resize_explicit_workspace(workspace, n, Solar, true);
            }
            Wide d_albedo = Wide::Zero();
            Wide d_thermal_surface = Wide::Zero();
            const Wide radiance = explicit_plane_view<Solar>(
                columns[0], views[view], albedo, thermal_surface, jacobians,
                d_albedo, d_thermal_surface, workspace);
            for (int lane = 0; lane < valid_lanes; ++lane) {
                worker.radiance[view * block_capacity + packet_base + lane] =
                    radiance[lane];
            }
            if (!jacobians) {
                continue;
            }

            const auto [bvp_albedo, bvp_thermal] =
                reverse_explicit_bvp<Solar>(columns[0], albedo, workspace);
            d_albedo += bvp_albedo;
            d_thermal_surface += bvp_thermal;
            reverse_explicit_layers<Solar>(columns[0], workspace);
            map_explicit_to_levels<Solar>(columns[0], workspace);

            for (int level = 0; level < nlevels; ++level) {
                const std::size_t offset =
                    (static_cast<std::size_t>(view) * nlevels + level) *
                        block_capacity +
                    packet_base;
                for (int lane = 0; lane < valid_lanes; ++lane) {
                    worker.extinction[offset + lane] =
                        workspace.d_level_extinction[level][lane];
                    worker.ssa[offset + lane] =
                        workspace.d_level_ssa[level][lane];
                    worker.b1[offset + lane] =
                        workspace.d_level_b1[level][lane];
                    if constexpr (!Solar) {
                        worker.emission[offset + lane] =
                            workspace.d_level_emission[level][lane];
                    }
                }
            }
            const std::size_t surface_offset =
                static_cast<std::size_t>(view) * block_capacity + packet_base;
            for (int lane = 0; lane < valid_lanes; ++lane) {
                worker.surface_albedo[surface_offset + lane] = d_albedo[lane];
                if constexpr (!Solar) {
                    worker.surface_emission[surface_offset + lane] =
                        d_thermal_surface[lane];
                }
            }
        }
        return true;
    }

    template <bool Solar>
    bool calculate_explicit_spherical_packet(
        const sasktran2::WavelengthBlock<>& block, int packet_base,
        Worker& worker) {
        if (!spherical || spherical->columns.empty()) {
            return false;
        }
        const int nlevels = geometry.size();
        const int n = nlevels - 1;
        const int nviews = static_cast<int>(spherical->ground_hit.size());
        const int valid_lanes = std::min(LANES, block.count - packet_base);
        const bool jacobians = atmosphere->num_deriv() > 0;
        auto& workspace = worker.explicit_spherical;
        workspace.columns.resize(spherical->columns.size());

        auto& first = workspace.columns.front();
        resize_explicit_workspace(first, n, Solar, false);
        for (int level = 0; level < nlevels; ++level) {
            first.level_extinction[level] =
                load_extinction(level, block, packet_base);
            first.level_ssa[level] = load_ssa(level, block, packet_base);
            first.level_b1[level] = load_b1(level, block, packet_base);
            if constexpr (!Solar) {
                first.level_emission[level] =
                    load_emission(level, block, packet_base);
            }
        }
        const Wide albedo = load_albedo(block, packet_base);
        Wide irradiance = Wide::Zero();
        Wide thermal_surface = Wide::Zero();
        if constexpr (Solar) {
            irradiance = load_spectral(atmosphere->storage().solar_irradiance,
                                       block, packet_base);
        } else {
            thermal_surface = load_spectral(atmosphere->surface().emission(),
                                            block, packet_base);
        }
        prepare_explicit_column<Solar>(spherical->columns.front(), irradiance,
                                       first);
        forward_explicit_layers<Solar>(spherical->columns.front(), albedo,
                                       thermal_surface, first);

        for (std::size_t index = 1; index < spherical->columns.size();
             ++index) {
            auto& column = workspace.columns[index];
            const auto& column_geometry = spherical->columns[index];
            resize_explicit_workspace(column, n, Solar, false);
            column.level_extinction = first.level_extinction;
            column.level_ssa = first.level_ssa;
            column.level_b1 = first.level_b1;
            column.od = first.od;
            column.ssa = first.ssa;
            column.b1 = first.b1;
            if constexpr (Solar) {
                column.attenuation[0] = Wide::Zero();
                column.transmission[0] = first.transmission[0];
                for (int boundary = 0; boundary < n; ++boundary) {
                    Wide slant = Wide::Zero();
                    for (int layer = 0; layer < n; ++layer) {
                        const double factor =
                            column_geometry.chapman(boundary, layer);
                        if (factor != 0.0) {
                            slant += column.od[layer] * factor;
                        }
                    }
                    column.attenuation[boundary + 1] = -slant;
                    column.transmission[boundary + 1] =
                        (-slant).exp() * column.transmission[0];
                }
                for (int layer = 0; layer < n; ++layer) {
                    column.secant[layer] =
                        positive_ratio(column.attenuation[layer] -
                                           column.attenuation[layer + 1],
                                       column.od[layer]);
                }
            } else {
                column.level_emission = first.level_emission;
                column.thermal_b0 = first.thermal_b0;
                column.thermal_b1 = first.thermal_b1;
            }
            forward_explicit_layers<Solar>(column_geometry, albedo,
                                           thermal_surface, column);
        }

        if (jacobians) {
            for (std::size_t index = 0; index < spherical->columns.size();
                 ++index) {
                if (!explicit_column_nonresonant<Solar>(
                        spherical->columns[index], workspace.columns[index])) {
                    return false;
                }
            }
        }

        const auto sza_interpolation = [this](double cosine) {
            const auto& grid = spherical->sza_grid;
            if (grid.size() == 1 || cosine <= grid.front()) {
                return std::tuple<std::size_t, std::size_t, double>{0, 0, 0.0};
            }
            const std::size_t last = grid.size() - 1;
            if (cosine >= grid.back()) {
                return std::tuple<std::size_t, std::size_t, double>{last, last,
                                                                    0.0};
            }
            const auto iterator =
                std::lower_bound(grid.begin(), grid.end(), cosine);
            const std::size_t upper =
                static_cast<std::size_t>(iterator - grid.begin());
            const std::size_t lower = upper - 1;
            const double weight =
                (cosine - grid[lower]) / (grid[upper] - grid[lower]);
            return std::tuple<std::size_t, std::size_t, double>{lower, upper,
                                                                weight};
        };

        const auto interpolated_local_source =
            [this, &workspace,
             &sza_interpolation](std::size_t segment) -> Wide {
            const auto [lower, upper, upper_weight] =
                sza_interpolation(spherical->segment_cos_sza[segment]);
            const auto evaluate = [this, &workspace,
                                   segment](std::size_t column) {
                return explicit_local_source<Solar>(
                    spherical->columns[column],
                    static_cast<int>(spherical->segment_layers[segment]),
                    spherical->segment_fractions[segment],
                    spherical->segment_cosines[segment],
                    spherical->segment_relative_azimuths[segment],
                    workspace.columns[column]);
            };
            const Wide lower_source = evaluate(lower);
            if (lower == upper) {
                return lower_source;
            }
            return lower_source * (1.0 - upper_weight) +
                   evaluate(upper) * upper_weight;
        };

        const auto interpolated_surface_source =
            [this, &workspace, &sza_interpolation, &albedo,
             &thermal_surface](int view) -> Wide {
            const auto [lower, upper, upper_weight] =
                sza_interpolation(spherical->ground_cos_sza[view]);
            const Wide lower_source = explicit_surface_source<Solar>(
                spherical->columns[lower], albedo, thermal_surface,
                workspace.columns[lower]);
            if (lower == upper) {
                return lower_source;
            }
            return lower_source * (1.0 - upper_weight) +
                   explicit_surface_source<Solar>(spherical->columns[upper],
                                                  albedo, thermal_surface,
                                                  workspace.columns[upper]) *
                       upper_weight;
        };

        const auto forward_ray =
            [this, &workspace, &interpolated_local_source,
             &interpolated_surface_source](int view) -> Wide {
            const std::size_t start = spherical->ray_offsets[view];
            const std::size_t end = spherical->ray_offsets[view + 1];
            const std::size_t segments = end - start;
            workspace.path_radiance.resize(segments + 1);
            workspace.path_transmission.resize(segments);
            workspace.path_source.resize(segments);
            workspace.path_radiance[0] = spherical->ground_hit[view]
                                             ? interpolated_surface_source(view)
                                             : Wide::Zero();
            for (std::size_t local = 0; local < segments; ++local) {
                const std::size_t segment = start + local;
                Wide optical_depth = Wide::Zero();
                const std::size_t od_start = spherical->od_offsets[segment];
                const std::size_t od_end = spherical->od_offsets[segment + 1];
                for (std::size_t stencil = od_start; stencil < od_end;
                     ++stencil) {
                    optical_depth +=
                        workspace.columns.front()
                            .level_extinction[spherical->od_indices[stencil]] *
                        spherical->od_weights[stencil];
                }
                const Wide transmission = (-optical_depth).exp();
                const Wide source = interpolated_local_source(segment);
                workspace.path_transmission[local] = transmission;
                workspace.path_source[local] = source;
                workspace.path_radiance[local + 1] =
                    workspace.path_radiance[local] * transmission +
                    source * (1.0 - transmission);
            }
            return workspace.path_radiance[segments];
        };

        for (int view = 0; view < nviews; ++view) {
            if (jacobians) {
                for (auto& column : workspace.columns) {
                    resize_explicit_workspace(column, n, Solar, true);
                }
                resize_wide(workspace.d_level_extinction, nlevels, true);
                resize_wide(workspace.d_level_ssa, nlevels, true);
                resize_wide(workspace.d_level_b1, nlevels, true);
                if constexpr (!Solar) {
                    resize_wide(workspace.d_level_emission, nlevels, true);
                }
            }

            const Wide radiance = forward_ray(view);
            for (int lane = 0; lane < valid_lanes; ++lane) {
                worker.radiance[view * block_capacity + packet_base + lane] =
                    radiance[lane];
            }
            if (!jacobians) {
                continue;
            }

            const std::size_t start = spherical->ray_offsets[view];
            const std::size_t end = spherical->ray_offsets[view + 1];
            Wide d_radiance = Wide::Ones();
            for (std::size_t local = end - start; local-- > 0;) {
                const std::size_t segment = start + local;
                const Wide& transmission = workspace.path_transmission[local];
                const Wide& source = workspace.path_source[local];
                const Wide& incoming = workspace.path_radiance[local];
                const Wide d_transmission = d_radiance * (incoming - source);
                const Wide d_source = d_radiance * (1.0 - transmission);
                const Wide d_optical_depth = -d_transmission * transmission;
                const std::size_t od_start = spherical->od_offsets[segment];
                const std::size_t od_end = spherical->od_offsets[segment + 1];
                for (std::size_t stencil = od_start; stencil < od_end;
                     ++stencil) {
                    workspace
                        .d_level_extinction[spherical->od_indices[stencil]] +=
                        d_optical_depth * spherical->od_weights[stencil];
                }

                const auto [lower, upper, upper_weight] =
                    sza_interpolation(spherical->segment_cos_sza[segment]);
                reverse_explicit_local_source<Solar>(
                    spherical->columns[lower],
                    static_cast<int>(spherical->segment_layers[segment]),
                    spherical->segment_fractions[segment],
                    spherical->segment_cosines[segment],
                    spherical->segment_relative_azimuths[segment],
                    d_source * (1.0 - upper_weight), workspace.columns[lower]);
                if (lower != upper) {
                    reverse_explicit_local_source<Solar>(
                        spherical->columns[upper],
                        static_cast<int>(spherical->segment_layers[segment]),
                        spherical->segment_fractions[segment],
                        spherical->segment_cosines[segment],
                        spherical->segment_relative_azimuths[segment],
                        d_source * upper_weight, workspace.columns[upper]);
                }
                d_radiance *= transmission;
            }

            Wide d_albedo = Wide::Zero();
            Wide d_thermal_surface = Wide::Zero();
            if (spherical->ground_hit[view]) {
                const auto [lower, upper, upper_weight] =
                    sza_interpolation(spherical->ground_cos_sza[view]);
                auto ground = reverse_explicit_surface_source<Solar>(
                    spherical->columns[lower], albedo,
                    d_radiance * (1.0 - upper_weight),
                    workspace.columns[lower]);
                d_albedo += ground.first;
                d_thermal_surface += ground.second;
                if (lower != upper) {
                    ground = reverse_explicit_surface_source<Solar>(
                        spherical->columns[upper], albedo,
                        d_radiance * upper_weight, workspace.columns[upper]);
                    d_albedo += ground.first;
                    d_thermal_surface += ground.second;
                }
            }

            for (std::size_t index = 0; index < spherical->columns.size();
                 ++index) {
                auto& column = workspace.columns[index];
                const auto boundary = reverse_explicit_bvp<Solar>(
                    spherical->columns[index], albedo, column);
                d_albedo += boundary.first;
                d_thermal_surface += boundary.second;
                reverse_explicit_layers<Solar>(spherical->columns[index],
                                               column);
                map_explicit_to_levels<Solar>(spherical->columns[index],
                                              column);
                for (int level = 0; level < nlevels; ++level) {
                    workspace.d_level_extinction[level] +=
                        column.d_level_extinction[level];
                    workspace.d_level_ssa[level] += column.d_level_ssa[level];
                    workspace.d_level_b1[level] += column.d_level_b1[level];
                    if constexpr (!Solar) {
                        workspace.d_level_emission[level] +=
                            column.d_level_emission[level];
                    }
                }
            }

            for (int level = 0; level < nlevels; ++level) {
                const std::size_t offset =
                    (static_cast<std::size_t>(view) * nlevels + level) *
                        block_capacity +
                    packet_base;
                for (int lane = 0; lane < valid_lanes; ++lane) {
                    worker.extinction[offset + lane] =
                        workspace.d_level_extinction[level][lane];
                    worker.ssa[offset + lane] =
                        workspace.d_level_ssa[level][lane];
                    worker.b1[offset + lane] =
                        workspace.d_level_b1[level][lane];
                    if constexpr (!Solar) {
                        worker.emission[offset + lane] =
                            workspace.d_level_emission[level][lane];
                    }
                }
            }
            const std::size_t surface_offset =
                static_cast<std::size_t>(view) * block_capacity + packet_base;
            for (int lane = 0; lane < valid_lanes; ++lane) {
                worker.surface_albedo[surface_offset + lane] = d_albedo[lane];
                if constexpr (!Solar) {
                    worker.surface_emission[surface_offset + lane] =
                        d_thermal_surface[lane];
                }
            }
        }
        return true;
    }

    void calculate_packet(const sasktran2::WavelengthBlock<>& block,
                          int packet_base, Worker& worker) {
        constexpr bool Solar = sasktran2::twostream::has_solar<SOURCE_TYPE>();
        const int nlevels = geometry.size();
        const int n = nlevels - 1;
        const int nviews = static_cast<int>(
            spherical ? spherical->ground_hit.size() : views.size());
        if (spherical && calculate_explicit_spherical_packet<Solar>(
                             block, packet_base, worker)) {
            return;
        }
        if (!spherical &&
            calculate_explicit_packet<Solar>(block, packet_base, worker)) {
            return;
        }
        auto& workspace = worker.packet;
        const std::size_t estimated_nodes =
            static_cast<std::size_t>(n) *
                (550 + 350 * std::max<std::size_t>(1, columns.size())) +
            static_cast<std::size_t>(nviews) * 300;
        const bool jacobians = atmosphere->num_deriv() > 0;
        workspace.tape.reset(estimated_nodes, jacobians);
        const Var zero = workspace.tape.constant(0.0);
        resize_vars(workspace.extinction, nlevels, zero);
        resize_vars(workspace.level_ssa, nlevels, zero);
        resize_vars(workspace.level_b1, nlevels, zero);
        resize_vars(workspace.emission, nlevels, zero);
        for (int level = 0; level < nlevels; ++level) {
            workspace.extinction[level] =
                workspace.tape.leaf(load_extinction(level, block, packet_base));
            workspace.level_ssa[level] =
                workspace.tape.leaf(load_ssa(level, block, packet_base));
            workspace.level_b1[level] =
                workspace.tape.leaf(load_b1(level, block, packet_base));
            if constexpr (!Solar) {
                workspace.emission[level] = workspace.tape.leaf(
                    load_emission(level, block, packet_base));
            }
        }
        const Var albedo = workspace.tape.leaf(load_albedo(block, packet_base));
        Var irradiance = zero;
        Var thermal_surface = zero;
        if constexpr (Solar) {
            irradiance = workspace.tape.leaf(load_spectral(
                atmosphere->storage().solar_irradiance, block, packet_base));
        } else {
            thermal_surface = workspace.tape.leaf(load_spectral(
                atmosphere->surface().emission(), block, packet_base));
        }

        workspace.columns.resize(columns.size());
        for (std::size_t index = 0; index < columns.size(); ++index) {
            prepare_and_solve_column<Solar>(
                columns[index], workspace.extinction, workspace.level_ssa,
                workspace.level_b1, workspace.emission, irradiance, albedo,
                thermal_surface, workspace.columns[index], zero);
        }

        workspace.outputs.clear();
        workspace.outputs.reserve(nviews);
        for (int view = 0; view < nviews; ++view) {
            if (spherical) {
                workspace.outputs.push_back(spherical_radiance<Solar>(
                    *spherical, workspace.columns, workspace.extinction, view,
                    albedo, thermal_surface, zero));
            } else {
                workspace.outputs.push_back(plane_radiance<Solar>(
                    columns[0], workspace.columns[0], views[view], albedo,
                    thermal_surface, zero));
            }
        }

        const int valid_lanes = std::min(LANES, block.count - packet_base);
        for (int view = 0; view < nviews; ++view) {
            const Wide& value = workspace.outputs[view].value();
            for (int lane = 0; lane < valid_lanes; ++lane) {
                worker.radiance[view * block_capacity + packet_base + lane] =
                    value[lane];
            }
            if (!jacobians) {
                continue;
            }
            workspace.tape.reverse(workspace.outputs[view].index);
            for (int level = 0; level < nlevels; ++level) {
                const std::size_t offset =
                    (static_cast<std::size_t>(view) * nlevels + level) *
                        block_capacity +
                    packet_base;
                const Wide& dext =
                    workspace.tape.adjoints[workspace.extinction[level].index];
                const Wide& dssa =
                    workspace.tape.adjoints[workspace.level_ssa[level].index];
                const Wide& db1 =
                    workspace.tape.adjoints[workspace.level_b1[level].index];
                const Wide* demission = nullptr;
                if constexpr (!Solar) {
                    demission = &workspace.tape
                                     .adjoints[workspace.emission[level].index];
                }
                for (int lane = 0; lane < valid_lanes; ++lane) {
                    worker.extinction[offset + lane] = dext[lane];
                    worker.ssa[offset + lane] = dssa[lane];
                    worker.b1[offset + lane] = db1[lane];
                    if constexpr (!Solar) {
                        worker.emission[offset + lane] = (*demission)[lane];
                    }
                }
            }
            const std::size_t surface_offset =
                static_cast<std::size_t>(view) * block_capacity + packet_base;
            const Wide& dalbedo = workspace.tape.adjoints[albedo.index];
            const Wide* dsurface_emission = nullptr;
            if constexpr (!Solar) {
                dsurface_emission =
                    &workspace.tape.adjoints[thermal_surface.index];
            }
            for (int lane = 0; lane < valid_lanes; ++lane) {
                worker.surface_albedo[surface_offset + lane] = dalbedo[lane];
                if constexpr (!Solar) {
                    worker.surface_emission[surface_offset + lane] =
                        (*dsurface_emission)[lane];
                }
            }
        }
    }
};

template <sasktran2::twostream::SourceType SOURCE_TYPE>
CppTwoStreamSourceAdapter<SOURCE_TYPE>::CppTwoStreamSourceAdapter(
    const sasktran2::Geometry1D& geometry)
    : m_impl(std::make_unique<Impl>(geometry)) {}

template <sasktran2::twostream::SourceType SOURCE_TYPE>
CppTwoStreamSourceAdapter<SOURCE_TYPE>::~CppTwoStreamSourceAdapter() = default;

template <sasktran2::twostream::SourceType SOURCE_TYPE>
void CppTwoStreamSourceAdapter<SOURCE_TYPE>::initialize_config(
    const sasktran2::Config& config) {
    m_impl->config = &config;
}

namespace {
    template <bool Solar>
    std::pair<Var, Var> lpsum(int az, double phase_mu, double phase_sine,
                              const Var& ssa, const Var& b1, const Var& zero) {
        if (az == 0) {
            return {0.5 * ssa * (1.0 - b1 * phase_mu),
                    0.5 * ssa * (1.0 + b1 * phase_mu)};
        }
        if constexpr (Solar) {
            const Var value = ssa * b1 * phase_sine;
            return {value, value};
        }
        return {zero, zero};
    }

    template <bool Solar>
    Var surface_source(const ColumnGeometry& geometry,
                       const ColumnSolution& solution, const Var& albedo,
                       const Var& thermal_surface, const Var& zero) {
        const int last = static_cast<int>(geometry.layer_thickness.size()) - 1;
        const Var base = solution.particular[0].gpb[last] +
                         solution.bvp[0].rhs[2 * last] *
                             solution.homogeneous[0].xp[last] *
                             solution.homogeneous[0].omega[last] +
                         solution.bvp[0].rhs[2 * last + 1] *
                             solution.homogeneous[0].xm[last];
        const Var diffuse = base * (2.0 * geometry.quadrature_cosine) * albedo;
        if constexpr (Solar) {
            return diffuse + zero;
        }
        return diffuse + thermal_surface;
    }

    template <bool Solar>
    Var local_source(const ColumnGeometry& geometry,
                     const ColumnSolution& solution, int layer,
                     double fraction_from_top, double propagation_cosine,
                     double relative_azimuth, const Var& zero) {
        const double mu = geometry.quadrature_cosine;
        constexpr int naz = Solar ? 2 : 1;
        const double phase_mu = propagation_cosine * mu;
        const double phase_sine =
            0.25 * std::sqrt(std::max(
                       0.0, (1.0 - propagation_cosine * propagation_cosine) *
                                (1.0 - mu * mu)));
        const std::array<double, 2> azimuth_weight = {
            1.0, std::cos(relative_azimuth)};
        const double fraction_from_bottom = 1.0 - fraction_from_top;
        Var source = zero;

        for (int az = 0; az < naz; ++az) {
            const auto& h = solution.homogeneous[az];
            const auto& p = solution.particular[az];
            const auto [lp, lm] =
                lpsum<Solar>(az, phase_mu, phase_sine, solution.ssa[layer],
                             solution.b1[layer], zero);
            const Var yp = lp * h.xp[layer] + lm * h.xm[layer];
            const Var ym = lp * h.xm[layer] + lm * h.xp[layer];
            const Var top_exponential =
                (-h.k[layer] * solution.od[layer] * fraction_from_top).exp();
            const Var bottom_exponential =
                (-h.k[layer] * solution.od[layer] * fraction_from_bottom).exp();
            Var value =
                solution.bvp[az].rhs[2 * layer] * yp * top_exponential +
                solution.bvp[az].rhs[2 * layer + 1] * ym * bottom_exponential;

            if constexpr (Solar) {
                const Var plus = solution.transmission[layer] *
                                 stable_exp_difference(
                                     h.k[layer], solution.secant[layer],
                                     solution.od[layer] * fraction_from_top);
                const Var minus =
                    solution.transmission[layer] * fraction_from_bottom *
                    stable_exp_difference(solution.secant[layer] *
                                              fraction_from_top,
                                          solution.secant[layer] +
                                              h.k[layer] * fraction_from_bottom,
                                          solution.od[layer]);
                value =
                    value + p.ap[layer] * yp * plus + p.am[layer] * ym * minus;
            } else {
                const Var plus = solution.thermal_b0[layer] *
                                 stable_exp_difference(
                                     h.k[layer], solution.thermal_b1[layer],
                                     solution.od[layer] * fraction_from_top);
                const Var minus =
                    solution.thermal_b0[layer] * fraction_from_bottom *
                    stable_exp_difference(solution.thermal_b1[layer] *
                                              fraction_from_top,
                                          solution.thermal_b1[layer] +
                                              h.k[layer] * fraction_from_bottom,
                                          solution.od[layer]);
                value = value + p.at[layer] * (yp * plus + ym * minus);
            }
            source = source + azimuth_weight[az] * value;
        }

        if constexpr (!Solar) {
            const Var thermal = (-solution.thermal_b1[layer] *
                                 solution.od[layer] * fraction_from_top)
                                    .exp();
            source = source + solution.thermal_b0[layer] * thermal *
                                  (1.0 - solution.ssa[layer]);
        }
        return source;
    }

    template <bool Solar>
    Var plane_radiance(const ColumnGeometry& geometry,
                       const ColumnSolution& solution, const ViewGeometry& view,
                       const Var& albedo, const Var& thermal_surface,
                       const Var& zero) {
        const int n = static_cast<int>(geometry.layer_thickness.size());
        constexpr int naz = Solar ? 2 : 1;
        const double mu = geometry.quadrature_cosine;
        const double inverse_view = 1.0 / view.cosine;
        const double phase_mu = view.cosine * mu;
        const double phase_sine =
            0.25 * std::sqrt(std::max(0.0, (1.0 - view.cosine * view.cosine) *
                                               (1.0 - mu * mu)));
        const std::array<double, 2> azimuth_weight = {
            1.0, std::cos(view.relative_azimuth)};
        Var attenuation = zero + 1.0;
        Var integrated = zero;

        for (int layer = 0; layer < n; ++layer) {
            const Var beam = (-solution.od[layer] * inverse_view).exp();
            const Var rate =
                Solar ? solution.secant[layer] : solution.thermal_b1[layer];
            const Var source_integral =
                inverse_view * stable_exp_difference(zero, rate + inverse_view,
                                                     solution.od[layer]);
            Var source = zero;
            for (int az = 0; az < naz; ++az) {
                const auto& h = solution.homogeneous[az];
                const auto& p = solution.particular[az];
                const auto [lp, lm] =
                    lpsum<Solar>(az, phase_mu, phase_sine, solution.ssa[layer],
                                 solution.b1[layer], zero);
                const Var yp = lp * h.xp[layer] + lm * h.xm[layer];
                const Var ym = lp * h.xm[layer] + lm * h.xp[layer];
                const Var hm =
                    inverse_view * stable_exp_difference(h.k[layer],
                                                         zero + inverse_view,
                                                         solution.od[layer]);
                const Var hp =
                    inverse_view *
                    stable_exp_difference(zero, h.k[layer] + inverse_view,
                                          solution.od[layer]);
                const Var dp_ratio =
                    inverse_view * stable_integrated_exp_difference(
                                       rate + h.k[layer], rate + inverse_view,
                                       solution.od[layer]);
                const Var dm_ratio =
                    inverse_view * stable_integrated_exp_difference(
                                       h.k[layer] + inverse_view,
                                       rate + inverse_view, solution.od[layer]);
                Var particular = zero;
                if constexpr (Solar) {
                    const Var dp = solution.transmission[layer] * dp_ratio;
                    const Var dm = solution.transmission[layer] * dm_ratio;
                    particular = p.ap[layer] * yp * dm + p.am[layer] * ym * dp;
                } else {
                    const Var dp = solution.thermal_b0[layer] * dp_ratio;
                    const Var dm = solution.thermal_b0[layer] * dm_ratio;
                    particular = p.at[layer] * (yp * dm + ym * dp);
                }
                source = source +
                         azimuth_weight[az] *
                             (solution.bvp[az].rhs[2 * layer] * yp * hp +
                              solution.bvp[az].rhs[2 * layer + 1] * ym * hm +
                              particular);
            }
            if constexpr (!Solar) {
                source = source + solution.thermal_b0[layer] * source_integral *
                                      (1.0 - solution.ssa[layer]);
            }
            integrated = integrated + source * attenuation;
            attenuation = attenuation * beam;
        }
        return integrated +
               attenuation * surface_source<Solar>(geometry, solution, albedo,
                                                   thermal_surface, zero);
    }

    std::tuple<int, int, double>
    sza_interpolation(const std::vector<double>& grid, double cosine) {
        if (grid.size() == 1 || cosine <= grid.front()) {
            return {0, 0, 0.0};
        }
        const int last = static_cast<int>(grid.size()) - 1;
        if (cosine >= grid.back()) {
            return {last, last, 0.0};
        }
        const auto upper_it =
            std::lower_bound(grid.begin(), grid.end(), cosine);
        const int upper = static_cast<int>(upper_it - grid.begin());
        const int lower = upper - 1;
        const double weight =
            (cosine - grid[lower]) / (grid[upper] - grid[lower]);
        return {lower, upper, weight};
    }

    template <bool Solar>
    Var spherical_radiance(const SphericalGeometry& spherical,
                           const std::vector<ColumnSolution>& columns,
                           const std::vector<Var>& extinction, int view,
                           const Var& albedo, const Var& thermal_surface,
                           const Var& zero) {
        const std::size_t start = spherical.ray_offsets[view];
        const std::size_t end = spherical.ray_offsets[view + 1];
        Var radiance = zero;
        if (spherical.ground_hit[view] != 0) {
            const auto [lower, upper, weight] = sza_interpolation(
                spherical.sza_grid, spherical.ground_cos_sza[view]);
            radiance =
                surface_source<Solar>(spherical.columns[lower], columns[lower],
                                      albedo, thermal_surface, zero) *
                (1.0 - weight);
            if (lower != upper) {
                radiance =
                    radiance + surface_source<Solar>(spherical.columns[upper],
                                                     columns[upper], albedo,
                                                     thermal_surface, zero) *
                                   weight;
            }
        }

        for (std::size_t segment = start; segment < end; ++segment) {
            Var optical_depth = zero;
            for (std::size_t stencil = spherical.od_offsets[segment];
                 stencil < spherical.od_offsets[segment + 1]; ++stencil) {
                optical_depth =
                    optical_depth + extinction[spherical.od_indices[stencil]] *
                                        spherical.od_weights[stencil];
            }
            const Var transmission = (-optical_depth).exp();
            const auto [lower, upper, weight] = sza_interpolation(
                spherical.sza_grid, spherical.segment_cos_sza[segment]);
            Var source =
                local_source<Solar>(
                    spherical.columns[lower], columns[lower],
                    static_cast<int>(spherical.segment_layers[segment]),
                    spherical.segment_fractions[segment],
                    spherical.segment_cosines[segment],
                    spherical.segment_relative_azimuths[segment], zero) *
                (1.0 - weight);
            if (lower != upper) {
                source =
                    source +
                    local_source<Solar>(
                        spherical.columns[upper], columns[upper],
                        static_cast<int>(spherical.segment_layers[segment]),
                        spherical.segment_fractions[segment],
                        spherical.segment_cosines[segment],
                        spherical.segment_relative_azimuths[segment], zero) *
                        weight;
            }
            radiance = radiance * transmission + source * (1.0 - transmission);
        }
        return radiance;
    }
} // namespace

namespace {
    void resize_wide(std::vector<Wide>& values, int size, bool zero) {
        values.resize(size, Wide::Zero());
        if (zero) {
            std::fill(values.begin(), values.end(), Wide::Zero());
        }
    }

    void resize_explicit_workspace(ExplicitWorkspace& w, int n, bool solar,
                                   bool clear_adjoints) {
        const int nlevels = n + 1;
        for (auto* values : {&w.level_extinction, &w.level_ssa, &w.level_b1,
                             &w.level_emission}) {
            resize_wide(*values, nlevels, false);
        }
        for (auto* values : {&w.od, &w.ssa, &w.b1, &w.beam, &w.source}) {
            resize_wide(*values, n, false);
        }
        resize_wide(w.attenuation, nlevels, false);
        if (solar) {
            resize_wide(w.transmission, nlevels, false);
            resize_wide(w.secant, n, false);
        } else {
            resize_wide(w.thermal_b0, n, false);
            resize_wide(w.thermal_b1, n, false);
        }

        for (auto* values : {&w.d_od, &w.d_ssa, &w.d_b1}) {
            resize_wide(*values, n, clear_adjoints);
        }
        for (auto* values : {&w.d_level_extinction, &w.d_level_ssa,
                             &w.d_level_b1, &w.d_level_emission}) {
            resize_wide(*values, nlevels, clear_adjoints);
        }
        if (solar) {
            resize_wide(w.d_transmission, nlevels, clear_adjoints);
            resize_wide(w.d_secant, n, clear_adjoints);
        } else {
            resize_wide(w.d_thermal_b0, n, clear_adjoints);
            resize_wide(w.d_thermal_b1, n, clear_adjoints);
        }

        const int naz = solar ? 2 : 1;
        for (int az = 0; az < naz; ++az) {
            auto& h = w.homogeneous[az];
            for (auto* values :
                 {&h.d, &h.s, &h.k, &h.xp, &h.xm, &h.omega, &h.norm}) {
                resize_wide(*values, n, false);
            }
            auto& dh = w.d_homogeneous[az];
            for (auto* values : {&dh.k, &dh.xp, &dh.xm, &dh.omega}) {
                resize_wide(*values, n, clear_adjoints);
            }
            auto& p = w.particular[az];
            for (auto* values :
                 {&p.qp, &p.qm, &p.ap, &p.am, &p.at, &p.exponential, &p.cp,
                  &p.cm, &p.gpt, &p.gpb, &p.gmt, &p.gmb}) {
                resize_wide(*values, n, false);
            }
            auto& dp = w.d_particular[az];
            for (auto* values :
                 {&dp.ap, &dp.am, &dp.at, &dp.gpt, &dp.gpb, &dp.gmt, &dp.gmb}) {
                resize_wide(*values, n, clear_adjoints);
            }
            auto& q = w.bvp[az];
            for (auto* values :
                 {&q.e, &q.c, &q.d, &q.a, &q.b, &q.rhs, &q.inverse_mu, &q.alpha,
                  &q.beta, &q.gamma, &q.z}) {
                resize_wide(*values, 2 * n, false);
            }
            resize_wide(w.d_solution[az], 2 * n, clear_adjoints);
        }
    }

    bool column_source_nonresonant(const Wide& rate, const Wide& eigenvalue,
                                   const Wide& optical_depth) {
        for (int lane = 0; lane < LANES; ++lane) {
            const double thickness = optical_depth[lane];
            if (std::abs((rate[lane] - eigenvalue[lane]) * thickness) <=
                    1.0e-5 ||
                std::abs((rate[lane] + eigenvalue[lane]) * thickness) <=
                    1.0e-5) {
                return false;
            }
        }
        return true;
    }

    bool plane_source_nonresonant(const Wide& rate, const Wide& eigenvalue,
                                  const Wide& optical_depth,
                                  double inverse_view) {
        for (int lane = 0; lane < LANES; ++lane) {
            const double thickness = optical_depth[lane];
            const std::array<double, 5> differences = {
                rate[lane] + inverse_view, inverse_view - eigenvalue[lane],
                eigenvalue[lane] + inverse_view, rate[lane] + eigenvalue[lane],
                rate[lane] - eigenvalue[lane]};
            for (double difference : differences) {
                if (std::abs(difference * thickness) <= 1.0e-5) {
                    return false;
                }
            }
        }
        return true;
    }

    void pentadiagonal_solve(ExplicitBvp& q) {
        const int n = static_cast<int>(q.d.size());
        q.inverse_mu[0] = 1.0 / q.d[0];
        q.alpha[0] = q.a[0] * q.inverse_mu[0];
        q.beta[0] = q.b[0] * q.inverse_mu[0];
        q.z[0] = q.rhs[0] * q.inverse_mu[0];
        if (n > 1) {
            q.gamma[1] = q.c[1];
            q.inverse_mu[1] = 1.0 / (q.d[1] - q.alpha[0] * q.gamma[1]);
            q.alpha[1] = (q.a[1] - q.beta[0] * q.gamma[1]) * q.inverse_mu[1];
            q.beta[1] = q.b[1] * q.inverse_mu[1];
            q.z[1] = (q.rhs[1] - q.z[0] * q.gamma[1]) * q.inverse_mu[1];
        }
        for (int index = 2; index < n; ++index) {
            q.gamma[index] = q.c[index] - q.alpha[index - 2] * q.e[index];
            q.inverse_mu[index] =
                1.0 / (q.d[index] - q.beta[index - 2] * q.e[index] -
                       q.alpha[index - 1] * q.gamma[index]);
            if (index + 1 < n) {
                q.alpha[index] =
                    (q.a[index] - q.beta[index - 1] * q.gamma[index]) *
                    q.inverse_mu[index];
            }
            if (index + 2 < n) {
                q.beta[index] = q.b[index] * q.inverse_mu[index];
            }
            q.z[index] = (q.rhs[index] - q.z[index - 2] * q.e[index] -
                          q.z[index - 1] * q.gamma[index]) *
                         q.inverse_mu[index];
        }
        q.rhs[n - 1] = q.z[n - 1];
        if (n > 1) {
            q.rhs[n - 2] = q.z[n - 2] - q.alpha[n - 2] * q.rhs[n - 1];
        }
        for (int index = n - 3; index >= 0; --index) {
            q.rhs[index] = q.z[index] - q.alpha[index] * q.rhs[index + 1] -
                           q.beta[index] * q.rhs[index + 2];
        }
    }

    void pentadiagonal_transpose_solve(const ExplicitBvp& q,
                                       std::vector<Wide>& rhs) {
        const int n = static_cast<int>(rhs.size());
        if (n > 1) {
            rhs[1] -= q.alpha[0] * rhs[0];
        }
        for (int index = 2; index < n; ++index) {
            rhs[index] -= q.alpha[index - 1] * rhs[index - 1] +
                          q.beta[index - 2] * rhs[index - 2];
        }
        rhs[n - 1] *= q.inverse_mu[n - 1];
        if (n > 1) {
            rhs[n - 2] = (rhs[n - 2] - q.gamma[n - 1] * rhs[n - 1]) *
                         q.inverse_mu[n - 2];
        }
        for (int index = n - 3; index >= 0; --index) {
            rhs[index] = (rhs[index] - q.gamma[index + 1] * rhs[index + 1] -
                          q.e[index + 2] * rhs[index + 2]) *
                         q.inverse_mu[index];
        }
    }

    template <bool Solar>
    void prepare_explicit_column(const ColumnGeometry& geometry,
                                 const Wide& irradiance, ExplicitWorkspace& w) {
        const int n = static_cast<int>(geometry.layer_thickness.size());
        for (int layer = 0; layer < n; ++layer) {
            const Wide scattering_top =
                w.level_extinction[layer] * w.level_ssa[layer];
            const Wide scattering_bottom =
                w.level_extinction[layer + 1] * w.level_ssa[layer + 1];
            const Wide average_extinction =
                0.5 *
                (w.level_extinction[layer] + w.level_extinction[layer + 1]);
            const Wide average_scattering =
                0.5 * (scattering_top + scattering_bottom);
            w.od[layer] = average_extinction * geometry.layer_thickness[layer];
            w.ssa[layer] =
                positive_ratio(average_scattering, average_extinction);
            for (int lane = 0; lane < LANES; ++lane) {
                w.ssa[layer][lane] = std::min(w.ssa[layer][lane], 1.0 - 1.0e-9);
            }
            w.b1[layer] = positive_ratio(
                0.5 * (scattering_top * w.level_b1[layer] +
                       scattering_bottom * w.level_b1[layer + 1]),
                average_scattering);
            if constexpr (!Solar) {
                w.thermal_b0[layer] = w.level_emission[layer];
                w.thermal_b1[layer] = thermal_profile_slope(
                    w.level_emission[layer], w.level_emission[layer + 1],
                    w.od[layer]);
            }
        }
        if constexpr (Solar) {
            w.attenuation[0] = Wide::Zero();
            w.transmission[0] = irradiance;
            for (int boundary = 0; boundary < n; ++boundary) {
                Wide slant = Wide::Zero();
                for (int layer = 0; layer < n; ++layer) {
                    const double factor = geometry.chapman(boundary, layer);
                    if (factor != 0.0) {
                        slant += w.od[layer] * factor;
                    }
                }
                w.attenuation[boundary + 1] = -slant;
                w.transmission[boundary + 1] = (-slant).exp() * irradiance;
            }
            for (int layer = 0; layer < n; ++layer) {
                w.secant[layer] = positive_ratio(w.attenuation[layer] -
                                                     w.attenuation[layer + 1],
                                                 w.od[layer]);
            }
        }
    }

    template <bool Solar>
    void build_and_solve_explicit_bvp(const ColumnGeometry& geometry, int az,
                                      const Wide& albedo,
                                      const Wide& thermal_surface,
                                      ExplicitWorkspace& w) {
        const int n = static_cast<int>(geometry.layer_thickness.size());
        const int size = 2 * n;
        const double mu = geometry.quadrature_cosine;
        const auto& h = w.homogeneous[az];
        const auto& p = w.particular[az];
        auto& q = w.bvp[az];
        q.e[0] = Wide::Zero();
        q.c[0] = Wide::Zero();
        q.b[0] = Wide::Zero();
        q.rhs[0] = -p.gpt[0];
        for (int layer = 0; layer < n - 1; ++layer) {
            q.rhs[2 * layer + 1] = p.gmt[layer + 1] - p.gmb[layer];
            q.rhs[2 * layer + 2] = p.gpt[layer + 1] - p.gpb[layer];
        }
        const int last = size - 1;
        const double delta = az == 0 ? 1.0 : 0.0;
        Wide direct_boundary_source;
        if constexpr (Solar) {
            direct_boundary_source = delta * geometry.solar_cosine * albedo /
                                     EIGEN_PI * w.transmission[n];
        } else {
            direct_boundary_source = thermal_surface;
        }
        q.rhs[last] = direct_boundary_source -
                      (p.gmb[n - 1] - 2.0 * delta * mu * albedo * p.gpb[n - 1]);
        q.d[0] = h.xp[0];
        q.a[0] = h.xm[0] * h.omega[0];
        for (int layer = 0; layer < n - 1; ++layer) {
            const int row = 2 * layer;
            q.e[row + 1] = Wide::Zero();
            q.c[row + 1] = h.xm[layer] * h.omega[layer];
            q.d[row + 1] = h.xp[layer];
            q.a[row + 1] = -h.xm[layer + 1];
            q.b[row + 1] = -h.xp[layer + 1] * h.omega[layer + 1];
            q.e[row + 2] = h.xp[layer] * h.omega[layer];
            q.c[row + 2] = h.xm[layer];
            q.d[row + 2] = -h.xp[layer + 1];
            q.a[row + 2] = -h.xm[layer + 1] * h.omega[layer + 1];
            q.b[row + 2] = Wide::Zero();
        }
        q.e[last] = Wide::Zero();
        q.c[last] = (h.xm[n - 1] - 2.0 * mu * albedo * delta * h.xp[n - 1]) *
                    h.omega[n - 1];
        q.d[last] = h.xp[n - 1] - 2.0 * mu * albedo * delta * h.xm[n - 1];
        q.a[last] = Wide::Zero();
        q.b[last] = Wide::Zero();
        pentadiagonal_solve(q);
    }

    template <bool Solar>
    void
    forward_explicit_layers(const ColumnGeometry& geometry, const Wide& albedo,
                            const Wide& thermal_surface, ExplicitWorkspace& w) {
        const int n = static_cast<int>(geometry.layer_thickness.size());
        constexpr int naz = Solar ? 2 : 1;
        const double mu = geometry.quadrature_cosine;
        const double angular = std::sqrt(
            std::max(0.0, (1.0 - mu * mu) * (1.0 - geometry.solar_cosine *
                                                       geometry.solar_cosine)));
        for (int layer = 0; layer < n; ++layer) {
            const Wide rate = Solar ? w.secant[layer] : w.thermal_b1[layer];
            const Wide exponential = (-rate * w.od[layer]).exp();
            for (int az = 0; az < naz; ++az) {
                w.particular[az].exponential[layer] = exponential;
            }
        }
        for (int az = 0; az < naz; ++az) {
            auto& h = w.homogeneous[az];
            auto& p = w.particular[az];
            for (int layer = 0; layer < n; ++layer) {
                if (az == 0) {
                    h.d[layer] = w.ssa[layer] * w.b1[layer] * mu - 1.0 / mu;
                    h.s[layer] = (w.ssa[layer] - 1.0) / mu;
                } else {
                    h.d[layer] = Wide::Constant(-1.0 / mu);
                    h.s[layer] =
                        (w.ssa[layer] * w.b1[layer] * (1.0 - mu * mu) - 2.0) /
                        (2.0 * mu);
                }
                h.k[layer] = (h.s[layer] * h.d[layer]).sqrt();
                const Wide s_over_k = h.s[layer] / h.k[layer];
                h.xp[layer] = 0.5 * (1.0 - s_over_k);
                h.xm[layer] = 0.5 * (1.0 + s_over_k);
                h.omega[layer] = (-h.k[layer] * w.od[layer]).exp();
                h.norm[layer] =
                    mu * (h.xp[layer].square() - h.xm[layer].square());

                if constexpr (Solar) {
                    if (az == 0) {
                        p.qp[layer] =
                            w.ssa[layer] *
                            (1.0 + w.b1[layer] * geometry.solar_cosine * mu) /
                            FOUR_PI;
                        p.qm[layer] =
                            w.ssa[layer] *
                            (1.0 - w.b1[layer] * geometry.solar_cosine * mu) /
                            FOUR_PI;
                    } else {
                        p.qp[layer] =
                            w.ssa[layer] * w.b1[layer] * angular / FOUR_PI;
                        p.qm[layer] = p.qp[layer];
                    }
                    p.ap[layer] = (p.qp[layer] * h.xp[layer] +
                                   p.qm[layer] * h.xm[layer]) /
                                  h.norm[layer];
                    p.am[layer] = (p.qm[layer] * h.xp[layer] +
                                   p.qp[layer] * h.xm[layer]) /
                                  h.norm[layer];
                } else {
                    p.at[layer] = (1.0 - w.ssa[layer]) *
                                  (h.xp[layer] + h.xm[layer]) / h.norm[layer];
                }

                const Wide rate = Solar ? w.secant[layer] : w.thermal_b1[layer];
                const Wide amplitude =
                    Solar ? w.transmission[layer] : w.thermal_b0[layer];
                Wide cp_ratio;
                Wide cm_ratio;
                if (column_source_nonresonant(rate, h.k[layer], w.od[layer])) {
                    cp_ratio = (h.omega[layer] - p.exponential[layer]) /
                               (rate - h.k[layer]);
                    cm_ratio = (1.0 - h.omega[layer] * p.exponential[layer]) /
                               (rate + h.k[layer]);
                } else {
                    cp_ratio =
                        stable_exp_difference(h.k[layer], rate, w.od[layer]);
                    cm_ratio = stable_exp_difference(
                        Wide::Zero(), rate + h.k[layer], w.od[layer]);
                }
                p.cp[layer] = amplitude * cp_ratio;
                p.cm[layer] = amplitude * cm_ratio;
                if constexpr (Solar) {
                    p.gpt[layer] = p.am[layer] * p.cm[layer] * h.xm[layer];
                    p.gpb[layer] = p.ap[layer] * p.cp[layer] * h.xp[layer];
                    p.gmt[layer] = p.am[layer] * p.cm[layer] * h.xp[layer];
                    p.gmb[layer] = p.ap[layer] * p.cp[layer] * h.xm[layer];
                } else {
                    p.gpt[layer] = p.at[layer] * p.cm[layer] * h.xm[layer];
                    p.gpb[layer] = p.at[layer] * p.cp[layer] * h.xp[layer];
                    p.gmt[layer] = p.at[layer] * p.cm[layer] * h.xp[layer];
                    p.gmb[layer] = p.at[layer] * p.cp[layer] * h.xm[layer];
                }
            }
        }
        for (int az = 0; az < naz; ++az) {
            build_and_solve_explicit_bvp<Solar>(geometry, az, albedo,
                                                thermal_surface, w);
        }
    }

    template <bool Solar>
    bool explicit_column_nonresonant(const ColumnGeometry& geometry,
                                     const ExplicitWorkspace& w) {
        const int n = static_cast<int>(geometry.layer_thickness.size());
        constexpr int naz = Solar ? 2 : 1;
        for (int layer = 0; layer < n; ++layer) {
            const Wide& rate = Solar ? w.secant[layer] : w.thermal_b1[layer];
            for (int az = 0; az < naz; ++az) {
                if (!column_source_nonresonant(rate, w.homogeneous[az].k[layer],
                                               w.od[layer])) {
                    return false;
                }
            }
        }
        return true;
    }

    template <bool Solar>
    bool explicit_packet_nonresonant(const ColumnGeometry& geometry,
                                     const std::vector<ViewGeometry>& views,
                                     const ExplicitWorkspace& w) {
        if (!explicit_column_nonresonant<Solar>(geometry, w)) {
            return false;
        }
        const int n = static_cast<int>(geometry.layer_thickness.size());
        constexpr int naz = Solar ? 2 : 1;
        for (int layer = 0; layer < n; ++layer) {
            const Wide& rate = Solar ? w.secant[layer] : w.thermal_b1[layer];
            for (int az = 0; az < naz; ++az) {
                for (const auto& view : views) {
                    if (!plane_source_nonresonant(
                            rate, w.homogeneous[az].k[layer], w.od[layer],
                            1.0 / view.cosine)) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    template <bool Solar>
    std::pair<Wide, Wide> explicit_lpsum(int az, double phase_mu,
                                         double phase_sine, const Wide& ssa,
                                         const Wide& b1) {
        if (az == 0) {
            return {0.5 * ssa * (1.0 - b1 * phase_mu),
                    0.5 * ssa * (1.0 + b1 * phase_mu)};
        }
        if constexpr (Solar) {
            const Wide value = ssa * b1 * phase_sine;
            return {value, value};
        }
        return {Wide::Zero(), Wide::Zero()};
    }

    template <bool Solar>
    Wide
    explicit_local_source(const ColumnGeometry& geometry, int layer,
                          double fraction_from_top, double propagation_cosine,
                          double relative_azimuth, const ExplicitWorkspace& w) {
        constexpr int naz = Solar ? 2 : 1;
        const double mu = geometry.quadrature_cosine;
        const double phase_mu = propagation_cosine * mu;
        const double phase_sine =
            0.25 * std::sqrt(std::max(
                       0.0, (1.0 - propagation_cosine * propagation_cosine) *
                                (1.0 - mu * mu)));
        const std::array<double, 2> azimuth_weight = {
            1.0, std::cos(relative_azimuth)};
        const double fraction_from_bottom = 1.0 - fraction_from_top;
        Wide source = Wide::Zero();

        for (int az = 0; az < naz; ++az) {
            const auto& h = w.homogeneous[az];
            const auto& p = w.particular[az];
            const auto [lp, lm] = explicit_lpsum<Solar>(
                az, phase_mu, phase_sine, w.ssa[layer], w.b1[layer]);
            const Wide yp = lp * h.xp[layer] + lm * h.xm[layer];
            const Wide ym = lp * h.xm[layer] + lm * h.xp[layer];
            const Wide top_exponential =
                (-h.k[layer] * w.od[layer] * fraction_from_top).exp();
            const Wide bottom_exponential =
                (-h.k[layer] * w.od[layer] * fraction_from_bottom).exp();
            Wide value = w.bvp[az].rhs[2 * layer] * yp * top_exponential +
                         w.bvp[az].rhs[2 * layer + 1] * ym * bottom_exponential;

            if constexpr (Solar) {
                const Wide beam =
                    (-w.secant[layer] * w.od[layer] * fraction_from_top).exp();
                const Wide plus =
                    w.transmission[layer] *
                    exp_difference_ratio(top_exponential, beam, h.k[layer],
                                         w.secant[layer],
                                         w.od[layer] * fraction_from_top);
                const Wide minus =
                    w.transmission[layer] * fraction_from_bottom *
                    exp_difference_ratio(
                        beam, p.exponential[layer] * bottom_exponential,
                        w.secant[layer] * fraction_from_top,
                        w.secant[layer] + h.k[layer] * fraction_from_bottom,
                        w.od[layer]);
                value += p.ap[layer] * yp * plus + p.am[layer] * ym * minus;
            } else {
                const Wide thermal =
                    (-w.thermal_b1[layer] * w.od[layer] * fraction_from_top)
                        .exp();
                const Wide plus =
                    w.thermal_b0[layer] *
                    exp_difference_ratio(top_exponential, thermal, h.k[layer],
                                         w.thermal_b1[layer],
                                         w.od[layer] * fraction_from_top);
                const Wide minus =
                    w.thermal_b0[layer] * fraction_from_bottom *
                    exp_difference_ratio(
                        thermal, p.exponential[layer] * bottom_exponential,
                        w.thermal_b1[layer] * fraction_from_top,
                        w.thermal_b1[layer] + h.k[layer] * fraction_from_bottom,
                        w.od[layer]);
                value += p.at[layer] * (yp * plus + ym * minus);
            }
            source += azimuth_weight[az] * value;
        }

        if constexpr (!Solar) {
            const Wide thermal =
                (-w.thermal_b1[layer] * w.od[layer] * fraction_from_top).exp();
            source += w.thermal_b0[layer] * thermal * (1.0 - w.ssa[layer]);
        }
        return source;
    }

    template <bool Solar>
    void reverse_explicit_local_source(const ColumnGeometry& geometry,
                                       int layer, double fraction_from_top,
                                       double propagation_cosine,
                                       double relative_azimuth,
                                       const Wide& seed, ExplicitWorkspace& w) {
        constexpr int naz = Solar ? 2 : 1;
        const double mu = geometry.quadrature_cosine;
        const double phase_mu = propagation_cosine * mu;
        const double phase_sine =
            0.25 * std::sqrt(std::max(
                       0.0, (1.0 - propagation_cosine * propagation_cosine) *
                                (1.0 - mu * mu)));
        const std::array<double, 2> azimuth_weight = {
            1.0, std::cos(relative_azimuth)};
        const double fraction_from_bottom = 1.0 - fraction_from_top;

        for (int az = 0; az < naz; ++az) {
            const auto& h = w.homogeneous[az];
            const auto& p = w.particular[az];
            const auto [lp, lm] = explicit_lpsum<Solar>(
                az, phase_mu, phase_sine, w.ssa[layer], w.b1[layer]);
            const Wide yp = lp * h.xp[layer] + lm * h.xm[layer];
            const Wide ym = lp * h.xm[layer] + lm * h.xp[layer];
            const Wide top_exponential =
                (-h.k[layer] * w.od[layer] * fraction_from_top).exp();
            const Wide bottom_exponential =
                (-h.k[layer] * w.od[layer] * fraction_from_bottom).exp();
            const Wide source_scale = seed * azimuth_weight[az];
            const Wide& top_solution = w.bvp[az].rhs[2 * layer];
            const Wide& bottom_solution = w.bvp[az].rhs[2 * layer + 1];
            Wide d_yp = source_scale * top_solution * top_exponential;
            Wide d_ym = source_scale * bottom_solution * bottom_exponential;
            const Wide d_top_exponential = source_scale * top_solution * yp;
            const Wide d_bottom_exponential =
                source_scale * bottom_solution * ym;
            w.d_solution[az][2 * layer] += source_scale * yp * top_exponential;
            w.d_solution[az][2 * layer + 1] +=
                source_scale * ym * bottom_exponential;
            Wide d_k = Wide::Zero();

            if constexpr (Solar) {
                const Wide beam =
                    (-w.secant[layer] * w.od[layer] * fraction_from_top).exp();
                const Wide plus_ratio = exp_difference_ratio(
                    top_exponential, beam, h.k[layer], w.secant[layer],
                    w.od[layer] * fraction_from_top);
                const Wide minus_ratio =
                    fraction_from_bottom *
                    exp_difference_ratio(
                        beam, p.exponential[layer] * bottom_exponential,
                        w.secant[layer] * fraction_from_top,
                        w.secant[layer] + h.k[layer] * fraction_from_bottom,
                        w.od[layer]);
                const Wide plus = w.transmission[layer] * plus_ratio;
                const Wide minus = w.transmission[layer] * minus_ratio;

                w.d_particular[az].ap[layer] += source_scale * yp * plus;
                d_yp += source_scale * p.ap[layer] * plus;
                const Wide d_plus = source_scale * p.ap[layer] * yp;
                w.d_particular[az].am[layer] += source_scale * ym * minus;
                d_ym += source_scale * p.am[layer] * minus;
                const Wide d_minus = source_scale * p.am[layer] * ym;

                w.d_transmission[layer] +=
                    d_plus * plus_ratio + d_minus * minus_ratio;
                const auto [d_plus_k, d_plus_rate, d_plus_thickness] =
                    exp_difference_ratio_adjoint(
                        top_exponential, beam, h.k[layer], w.secant[layer],
                        w.od[layer] * fraction_from_top,
                        d_plus * w.transmission[layer]);
                const auto [d_minus_a, d_minus_b, d_minus_thickness] =
                    exp_difference_ratio_adjoint(
                        beam, p.exponential[layer] * bottom_exponential,
                        w.secant[layer] * fraction_from_top,
                        w.secant[layer] + h.k[layer] * fraction_from_bottom,
                        w.od[layer],
                        d_minus * w.transmission[layer] * fraction_from_bottom);
                d_k += d_plus_k + d_minus_b * fraction_from_bottom;
                w.d_secant[layer] +=
                    d_plus_rate + d_minus_a * fraction_from_top + d_minus_b;
                w.d_od[layer] +=
                    d_plus_thickness * fraction_from_top + d_minus_thickness;
            } else {
                const Wide thermal =
                    (-w.thermal_b1[layer] * w.od[layer] * fraction_from_top)
                        .exp();
                const Wide plus_ratio = exp_difference_ratio(
                    top_exponential, thermal, h.k[layer], w.thermal_b1[layer],
                    w.od[layer] * fraction_from_top);
                const Wide minus_ratio =
                    fraction_from_bottom *
                    exp_difference_ratio(
                        thermal, p.exponential[layer] * bottom_exponential,
                        w.thermal_b1[layer] * fraction_from_top,
                        w.thermal_b1[layer] + h.k[layer] * fraction_from_bottom,
                        w.od[layer]);
                const Wide plus = w.thermal_b0[layer] * plus_ratio;
                const Wide minus = w.thermal_b0[layer] * minus_ratio;

                w.d_particular[az].at[layer] +=
                    source_scale * (yp * plus + ym * minus);
                d_yp += source_scale * p.at[layer] * plus;
                d_ym += source_scale * p.at[layer] * minus;
                const Wide d_plus = source_scale * p.at[layer] * yp;
                const Wide d_minus = source_scale * p.at[layer] * ym;

                w.d_thermal_b0[layer] +=
                    d_plus * plus_ratio + d_minus * minus_ratio;
                const auto [d_plus_k, d_plus_rate, d_plus_thickness] =
                    exp_difference_ratio_adjoint(
                        top_exponential, thermal, h.k[layer],
                        w.thermal_b1[layer], w.od[layer] * fraction_from_top,
                        d_plus * w.thermal_b0[layer]);
                const auto [d_minus_a, d_minus_b, d_minus_thickness] =
                    exp_difference_ratio_adjoint(
                        thermal, p.exponential[layer] * bottom_exponential,
                        w.thermal_b1[layer] * fraction_from_top,
                        w.thermal_b1[layer] + h.k[layer] * fraction_from_bottom,
                        w.od[layer],
                        d_minus * w.thermal_b0[layer] * fraction_from_bottom);
                d_k += d_plus_k + d_minus_b * fraction_from_bottom;
                w.d_thermal_b1[layer] +=
                    d_plus_rate + d_minus_a * fraction_from_top + d_minus_b;
                w.d_od[layer] +=
                    d_plus_thickness * fraction_from_top + d_minus_thickness;
            }

            const Wide d_lp = d_yp * h.xp[layer] + d_ym * h.xm[layer];
            const Wide d_lm = d_yp * h.xm[layer] + d_ym * h.xp[layer];
            w.d_homogeneous[az].xp[layer] += d_yp * lp + d_ym * lm;
            w.d_homogeneous[az].xm[layer] += d_yp * lm + d_ym * lp;
            if (az == 0) {
                w.d_ssa[layer] += d_lp * 0.5 * (1.0 - w.b1[layer] * phase_mu) +
                                  d_lm * 0.5 * (1.0 + w.b1[layer] * phase_mu);
                w.d_b1[layer] += d_lp * (-0.5 * w.ssa[layer] * phase_mu) +
                                 d_lm * (0.5 * w.ssa[layer] * phase_mu);
            } else if constexpr (Solar) {
                w.d_ssa[layer] += (d_lp + d_lm) * w.b1[layer] * phase_sine;
                w.d_b1[layer] += (d_lp + d_lm) * w.ssa[layer] * phase_sine;
            }

            d_k -= d_top_exponential * top_exponential * w.od[layer] *
                   fraction_from_top;
            w.d_od[layer] -= d_top_exponential * top_exponential * h.k[layer] *
                             fraction_from_top;
            d_k -= d_bottom_exponential * bottom_exponential * w.od[layer] *
                   fraction_from_bottom;
            w.d_od[layer] -= d_bottom_exponential * bottom_exponential *
                             h.k[layer] * fraction_from_bottom;
            w.d_homogeneous[az].k[layer] += d_k;
        }

        if constexpr (!Solar) {
            const Wide thermal =
                (-w.thermal_b1[layer] * w.od[layer] * fraction_from_top).exp();
            w.d_thermal_b0[layer] += seed * thermal * (1.0 - w.ssa[layer]);
            const Wide d_thermal =
                seed * w.thermal_b0[layer] * (1.0 - w.ssa[layer]);
            w.d_ssa[layer] -= seed * w.thermal_b0[layer] * thermal;
            w.d_thermal_b1[layer] -=
                d_thermal * thermal * w.od[layer] * fraction_from_top;
            w.d_od[layer] -=
                d_thermal * thermal * w.thermal_b1[layer] * fraction_from_top;
        }
    }

    template <bool Solar>
    Wide explicit_surface_source(const ColumnGeometry& geometry,
                                 const Wide& albedo,
                                 const Wide& thermal_surface,
                                 const ExplicitWorkspace& w) {
        const int last = static_cast<int>(geometry.layer_thickness.size()) - 1;
        const Wide base =
            w.particular[0].gpb[last] +
            w.bvp[0].rhs[2 * last] * w.homogeneous[0].xp[last] *
                w.homogeneous[0].omega[last] +
            w.bvp[0].rhs[2 * last + 1] * w.homogeneous[0].xm[last];
        Wide result = base * (2.0 * geometry.quadrature_cosine) * albedo;
        if constexpr (!Solar) {
            result += thermal_surface;
        }
        return result;
    }

    template <bool Solar>
    std::pair<Wide, Wide>
    reverse_explicit_surface_source(const ColumnGeometry& geometry,
                                    const Wide& albedo, const Wide& seed,
                                    ExplicitWorkspace& w) {
        const int last = static_cast<int>(geometry.layer_thickness.size()) - 1;
        const Wide base =
            w.particular[0].gpb[last] +
            w.bvp[0].rhs[2 * last] * w.homogeneous[0].xp[last] *
                w.homogeneous[0].omega[last] +
            w.bvp[0].rhs[2 * last + 1] * w.homogeneous[0].xm[last];
        const double scale = 2.0 * geometry.quadrature_cosine;
        const Wide d_base = seed * scale * albedo;
        w.d_particular[0].gpb[last] += d_base;
        w.d_solution[0][2 * last] +=
            d_base * w.homogeneous[0].xp[last] * w.homogeneous[0].omega[last];
        w.d_homogeneous[0].xp[last] +=
            d_base * w.bvp[0].rhs[2 * last] * w.homogeneous[0].omega[last];
        w.d_homogeneous[0].omega[last] +=
            d_base * w.bvp[0].rhs[2 * last] * w.homogeneous[0].xp[last];
        w.d_solution[0][2 * last + 1] += d_base * w.homogeneous[0].xm[last];
        w.d_homogeneous[0].xm[last] += d_base * w.bvp[0].rhs[2 * last + 1];
        return {seed * scale * base, Solar ? Wide::Zero() : seed};
    }

    template <bool Solar>
    void reverse_explicit_view_azimuth_fast(
        int az, const ViewGeometry& view, double azimuth_weight,
        double phase_mu, double phase_sine, int layer, const Wide& beam,
        const Wide& exponential, const Wide& source_integral,
        const Wide& d_source, Wide& d_beam, Wide& d_exponential,
        Wide& d_source_integral, ExplicitWorkspace& w) {
        const double vmu = view.cosine;
        const auto& h = w.homogeneous[az];
        const auto& p = w.particular[az];
        const auto [lp, lm] = explicit_lpsum<Solar>(az, phase_mu, phase_sine,
                                                    w.ssa[layer], w.b1[layer]);
        const Wide yp = lp * h.xp[layer] + lm * h.xm[layer];
        const Wide ym = lp * h.xm[layer] + lm * h.xp[layer];
        const Wide hm_denominator = 1.0 - h.k[layer] * vmu;
        const Wide hp_denominator = 1.0 + h.k[layer] * vmu;
        const Wide inverse_hm = 1.0 / hm_denominator;
        const Wide inverse_hp = 1.0 / hp_denominator;
        const Wide hm = (h.omega[layer] - beam) * inverse_hm;
        const Wide hp = (1.0 - h.omega[layer] * beam) * inverse_hp;
        const Wide source_scale = d_source * azimuth_weight;
        Wide d_yp = source_scale * w.bvp[az].rhs[2 * layer] * hp;
        Wide d_ym = source_scale * w.bvp[az].rhs[2 * layer + 1] * hm;
        Wide d_hp = source_scale * w.bvp[az].rhs[2 * layer] * yp;
        Wide d_hm = source_scale * w.bvp[az].rhs[2 * layer + 1] * ym;
        w.d_solution[az][2 * layer] += source_scale * yp * hp;
        w.d_solution[az][2 * layer + 1] += source_scale * ym * hm;

        const Wide dp_denominator =
            (Solar ? w.secant[layer] : w.thermal_b1[layer]) + h.k[layer];
        const Wide dm_denominator =
            (Solar ? w.secant[layer] : w.thermal_b1[layer]) - h.k[layer];
        const Wide inverse_dp = 1.0 / dp_denominator;
        const Wide inverse_dm = 1.0 / dm_denominator;
        const Wide amplitude =
            Solar ? w.transmission[layer] : w.thermal_b0[layer];
        const Wide dp =
            amplitude * (source_integral - exponential * hm) * inverse_dp;
        const Wide dm = amplitude * (hp - source_integral) * inverse_dm;
        Wide d_dp;
        Wide d_dm;
        if constexpr (Solar) {
            w.d_particular[az].ap[layer] += source_scale * yp * dm;
            d_yp += source_scale * p.ap[layer] * dm;
            d_dm = source_scale * p.ap[layer] * yp;
            w.d_particular[az].am[layer] += source_scale * ym * dp;
            d_ym += source_scale * p.am[layer] * dp;
            d_dp = source_scale * p.am[layer] * ym;
        } else {
            w.d_particular[az].at[layer] += source_scale * (yp * dm + ym * dp);
            d_yp += source_scale * p.at[layer] * dm;
            d_ym += source_scale * p.at[layer] * dp;
            d_dm = source_scale * p.at[layer] * yp;
            d_dp = source_scale * p.at[layer] * ym;
        }

        const Wide dp_numerator = source_integral - exponential * hm;
        const Wide d_dp_numerator = d_dp * amplitude * inverse_dp;
        const Wide d_dp_denominator = -d_dp * dp * inverse_dp;
        const Wide dm_numerator = hp - source_integral;
        const Wide d_dm_numerator = d_dm * amplitude * inverse_dm;
        const Wide d_dm_denominator = -d_dm * dm * inverse_dm;
        if constexpr (Solar) {
            w.d_transmission[layer] += d_dp * dp_numerator * inverse_dp +
                                       d_dm * dm_numerator * inverse_dm;
            w.d_secant[layer] += d_dp_denominator + d_dm_denominator;
        } else {
            w.d_thermal_b0[layer] += d_dp * dp_numerator * inverse_dp +
                                     d_dm * dm_numerator * inverse_dm;
            w.d_thermal_b1[layer] += d_dp_denominator + d_dm_denominator;
        }
        d_source_integral += d_dp_numerator - d_dm_numerator;
        d_exponential -= d_dp_numerator * hm;
        d_hm -= d_dp_numerator * exponential;
        d_hp += d_dm_numerator;
        w.d_homogeneous[az].k[layer] += d_dp_denominator - d_dm_denominator;

        const Wide d_hp_numerator = d_hp * inverse_hp;
        const Wide d_hp_denominator = -d_hp * hp * inverse_hp;
        w.d_homogeneous[az].omega[layer] -= d_hp_numerator * beam;
        d_beam -= d_hp_numerator * h.omega[layer];
        w.d_homogeneous[az].k[layer] += d_hp_denominator * vmu;

        const Wide d_hm_numerator = d_hm * inverse_hm;
        const Wide d_hm_denominator = -d_hm * hm * inverse_hm;
        w.d_homogeneous[az].omega[layer] += d_hm_numerator;
        d_beam -= d_hm_numerator;
        w.d_homogeneous[az].k[layer] -= d_hm_denominator * vmu;

        const Wide d_lp = d_yp * h.xp[layer] + d_ym * h.xm[layer];
        const Wide d_lm = d_yp * h.xm[layer] + d_ym * h.xp[layer];
        w.d_homogeneous[az].xp[layer] += d_yp * lp + d_ym * lm;
        w.d_homogeneous[az].xm[layer] += d_yp * lm + d_ym * lp;
        if (az == 0) {
            w.d_ssa[layer] += d_lp * 0.5 * (1.0 - w.b1[layer] * phase_mu) +
                              d_lm * 0.5 * (1.0 + w.b1[layer] * phase_mu);
            w.d_b1[layer] += d_lp * (-0.5 * w.ssa[layer] * phase_mu) +
                             d_lm * (0.5 * w.ssa[layer] * phase_mu);
        } else if constexpr (Solar) {
            w.d_ssa[layer] += (d_lp + d_lm) * w.b1[layer] * phase_sine;
            w.d_b1[layer] += (d_lp + d_lm) * w.ssa[layer] * phase_sine;
        }
    }

    template <bool Solar>
    Wide explicit_plane_view(const ColumnGeometry& geometry,
                             const ViewGeometry& view, const Wide& albedo,
                             const Wide& thermal_surface, bool reverse,
                             Wide& d_albedo, Wide& d_thermal_surface,
                             ExplicitWorkspace& w) {
        const int n = static_cast<int>(geometry.layer_thickness.size());
        constexpr int naz = Solar ? 2 : 1;
        const double mu = geometry.quadrature_cosine;
        const double inverse_view = 1.0 / view.cosine;
        const double phase_mu = view.cosine * mu;
        const double phase_sine =
            0.25 * std::sqrt(std::max(0.0, (1.0 - view.cosine * view.cosine) *
                                               (1.0 - mu * mu)));
        const std::array<double, 2> azimuth_weight = {
            1.0, std::cos(view.relative_azimuth)};
        bool fast = true;
        for (int layer = 0; layer < n; ++layer) {
            const Wide& rate = Solar ? w.secant[layer] : w.thermal_b1[layer];
            for (int az = 0; az < naz; ++az) {
                fast = fast && plane_source_nonresonant(
                                   rate, w.homogeneous[az].k[layer],
                                   w.od[layer], inverse_view);
            }
        }
        if (reverse && !fast) {
            throw std::logic_error(
                "explicit two-stream reverse called for resonant packet");
        }

        w.attenuation[0] = Wide::Ones();
        Wide integrated = Wide::Zero();
        for (int layer = 0; layer < n; ++layer) {
            const Wide beam = (-w.od[layer] * inverse_view).exp();
            w.beam[layer] = beam;
            const Wide& exponential = w.particular[0].exponential[layer];
            const Wide& rate = Solar ? w.secant[layer] : w.thermal_b1[layer];
            Wide source_integral;
            if (fast) {
                source_integral =
                    (1.0 - exponential * beam) / (1.0 + rate * view.cosine);
            } else {
                source_integral =
                    inverse_view *
                    stable_exp_difference(Wide::Zero(),
                                          rate + Wide::Constant(inverse_view),
                                          w.od[layer]);
            }
            w.source[layer] = Wide::Zero();
            for (int az = 0; az < naz; ++az) {
                const auto& h = w.homogeneous[az];
                const auto& p = w.particular[az];
                const auto [lp, lm] = explicit_lpsum<Solar>(
                    az, phase_mu, phase_sine, w.ssa[layer], w.b1[layer]);
                const Wide yp = lp * h.xp[layer] + lm * h.xm[layer];
                const Wide ym = lp * h.xm[layer] + lm * h.xp[layer];
                Wide hm;
                Wide hp;
                Wide dp_ratio;
                Wide dm_ratio;
                if (fast) {
                    hm = (h.omega[layer] - beam) /
                         (1.0 - h.k[layer] * view.cosine);
                    hp = (1.0 - h.omega[layer] * beam) /
                         (1.0 + h.k[layer] * view.cosine);
                    dp_ratio = (source_integral - exponential * hm) /
                               (rate + h.k[layer]);
                    dm_ratio = (hp - source_integral) / (rate - h.k[layer]);
                } else {
                    const Wide inv = Wide::Constant(inverse_view);
                    hm = inverse_view *
                         stable_exp_difference(h.k[layer], inv, w.od[layer]);
                    hp = inverse_view * stable_exp_difference(Wide::Zero(),
                                                              h.k[layer] + inv,
                                                              w.od[layer]);
                    dp_ratio = inverse_view *
                               stable_integrated_exp_difference(
                                   rate + h.k[layer], rate + inv, w.od[layer]);
                    dm_ratio = inverse_view *
                               stable_integrated_exp_difference(
                                   h.k[layer] + inv, rate + inv, w.od[layer]);
                }
                Wide particular;
                if constexpr (Solar) {
                    particular =
                        p.ap[layer] * yp * (w.transmission[layer] * dm_ratio) +
                        p.am[layer] * ym * (w.transmission[layer] * dp_ratio);
                } else {
                    particular =
                        p.at[layer] * (yp * (w.thermal_b0[layer] * dm_ratio) +
                                       ym * (w.thermal_b0[layer] * dp_ratio));
                }
                w.source[layer] +=
                    azimuth_weight[az] *
                    (w.bvp[az].rhs[2 * layer] * yp * hp +
                     w.bvp[az].rhs[2 * layer + 1] * ym * hm + particular);
            }
            if constexpr (!Solar) {
                w.source[layer] += w.thermal_b0[layer] * source_integral *
                                   (1.0 - w.ssa[layer]);
            }
            integrated += w.source[layer] * w.attenuation[layer];
            w.attenuation[layer + 1] = w.attenuation[layer] * beam;
        }

        const int last = n - 1;
        const Wide base_surface =
            w.particular[0].gpb[last] +
            w.bvp[0].rhs[2 * last] * w.homogeneous[0].xp[last] *
                w.homogeneous[0].omega[last] +
            w.bvp[0].rhs[2 * last + 1] * w.homogeneous[0].xm[last];
        const Wide surface = base_surface * (2.0 * mu) * albedo;
        const Wide output =
            integrated + w.attenuation[n] * (surface + thermal_surface);
        if (!reverse) {
            return output;
        }

        const Wide seed = Wide::Ones();
        const Wide d_surface = seed * w.attenuation[n];
        if constexpr (!Solar) {
            d_thermal_surface += d_surface;
        }
        const Wide d_base_surface = d_surface * (2.0 * mu) * albedo;
        d_albedo += d_surface * (2.0 * mu) * base_surface;
        w.d_particular[0].gpb[last] += d_base_surface;
        w.d_solution[0][2 * last] += d_base_surface *
                                     w.homogeneous[0].xp[last] *
                                     w.homogeneous[0].omega[last];
        w.d_homogeneous[0].xp[last] += d_base_surface * w.bvp[0].rhs[2 * last] *
                                       w.homogeneous[0].omega[last];
        w.d_homogeneous[0].omega[last] +=
            d_base_surface * w.bvp[0].rhs[2 * last] * w.homogeneous[0].xp[last];
        w.d_solution[0][2 * last + 1] +=
            d_base_surface * w.homogeneous[0].xm[last];
        w.d_homogeneous[0].xm[last] +=
            d_base_surface * w.bvp[0].rhs[2 * last + 1];

        Wide d_attenuation = seed * (surface + thermal_surface);
        for (int layer = n - 1; layer >= 0; --layer) {
            const Wide& beam = w.beam[layer];
            Wide d_beam = d_attenuation * w.attenuation[layer];
            d_attenuation = d_attenuation * beam + seed * w.source[layer];
            const Wide d_source = seed * w.attenuation[layer];
            const Wide& exponential = w.particular[0].exponential[layer];
            const Wide& rate = Solar ? w.secant[layer] : w.thermal_b1[layer];
            const Wide source_integral =
                (1.0 - exponential * beam) / (1.0 + rate * view.cosine);
            Wide d_source_integral = Wide::Zero();
            Wide d_exponential = Wide::Zero();
            for (int az = 0; az < naz; ++az) {
                reverse_explicit_view_azimuth_fast<Solar>(
                    az, view, azimuth_weight[az], phase_mu, phase_sine, layer,
                    beam, exponential, source_integral, d_source, d_beam,
                    d_exponential, d_source_integral, w);
            }
            if constexpr (!Solar) {
                w.d_thermal_b0[layer] +=
                    d_source * source_integral * (1.0 - w.ssa[layer]);
                d_source_integral +=
                    d_source * w.thermal_b0[layer] * (1.0 - w.ssa[layer]);
                w.d_ssa[layer] -=
                    d_source * w.thermal_b0[layer] * source_integral;
            }
            const Wide denominator = 1.0 + rate * view.cosine;
            const Wide inverse_denominator = 1.0 / denominator;
            const Wide d_numerator = d_source_integral * inverse_denominator;
            const Wide d_denominator =
                -d_source_integral * source_integral * inverse_denominator;
            d_exponential -= d_numerator * beam;
            d_beam -= d_numerator * exponential;
            if constexpr (Solar) {
                w.d_secant[layer] += d_denominator * view.cosine;
                w.d_secant[layer] -= d_exponential * exponential * w.od[layer];
                w.d_od[layer] -= d_exponential * exponential * w.secant[layer];
            } else {
                w.d_thermal_b1[layer] += d_denominator * view.cosine;
                w.d_thermal_b1[layer] -=
                    d_exponential * exponential * w.od[layer];
                w.d_od[layer] -=
                    d_exponential * exponential * w.thermal_b1[layer];
            }
            w.d_od[layer] -= d_beam * beam * inverse_view;
        }
        return output;
    }

    template <bool Solar>
    std::pair<Wide, Wide> reverse_explicit_bvp(const ColumnGeometry& geometry,
                                               const Wide& albedo,
                                               ExplicitWorkspace& w) {
        const int n = static_cast<int>(geometry.layer_thickness.size());
        const int size = 2 * n;
        constexpr int naz = Solar ? 2 : 1;
        const double mu = geometry.quadrature_cosine;
        Wide d_albedo = Wide::Zero();
        Wide d_thermal_surface = Wide::Zero();
        for (int az = 0; az < naz; ++az) {
            pentadiagonal_transpose_solve(w.bvp[az], w.d_solution[az]);
            const auto& lambda = w.d_solution[az];
            const auto& solution = w.bvp[az].rhs;
            const double delta = az == 0 ? 1.0 : 0.0;

            w.d_particular[az].gpt[0] -= lambda[0];
            for (int layer = 0; layer < n - 1; ++layer) {
                w.d_particular[az].gmt[layer + 1] += lambda[2 * layer + 1];
                w.d_particular[az].gmb[layer] -= lambda[2 * layer + 1];
                w.d_particular[az].gpt[layer + 1] += lambda[2 * layer + 2];
                w.d_particular[az].gpb[layer] -= lambda[2 * layer + 2];
            }
            const int last_row = size - 1;
            d_albedo += lambda[last_row] * 2.0 * delta * mu *
                        w.particular[az].gpb[n - 1];
            w.d_particular[az].gmb[n - 1] -= lambda[last_row];
            w.d_particular[az].gpb[n - 1] +=
                lambda[last_row] * 2.0 * delta * mu * albedo;
            if constexpr (Solar) {
                w.d_transmission[n] += lambda[last_row] * delta *
                                       geometry.solar_cosine * albedo /
                                       EIGEN_PI;
                d_albedo += lambda[last_row] * delta * geometry.solar_cosine /
                            EIGEN_PI * w.transmission[n];
            } else {
                d_thermal_surface += lambda[last_row];
            }

            const Wide g_d0 = -lambda[0] * solution[0];
            const Wide g_a0 = -lambda[0] * solution[1];
            w.d_homogeneous[az].xp[0] += g_d0;
            w.d_homogeneous[az].xm[0] += g_a0 * w.homogeneous[az].omega[0];
            w.d_homogeneous[az].omega[0] += g_a0 * w.homogeneous[az].xm[0];

            for (int layer = 0; layer < n - 1; ++layer) {
                int row = 2 * layer + 1;
                const Wide g_c = -lambda[row] * solution[row - 1];
                const Wide g_d = -lambda[row] * solution[row];
                const Wide g_a = -lambda[row] * solution[row + 1];
                const Wide g_b = -lambda[row] * solution[row + 2];
                w.d_homogeneous[az].xm[layer] +=
                    g_c * w.homogeneous[az].omega[layer];
                w.d_homogeneous[az].omega[layer] +=
                    g_c * w.homogeneous[az].xm[layer];
                w.d_homogeneous[az].xp[layer] += g_d;
                w.d_homogeneous[az].xm[layer + 1] -= g_a;
                w.d_homogeneous[az].xp[layer + 1] -=
                    g_b * w.homogeneous[az].omega[layer + 1];
                w.d_homogeneous[az].omega[layer + 1] -=
                    g_b * w.homogeneous[az].xp[layer + 1];

                row = 2 * layer + 2;
                const Wide g_e = -lambda[row] * solution[row - 2];
                const Wide g_c2 = -lambda[row] * solution[row - 1];
                const Wide g_d2 = -lambda[row] * solution[row];
                Wide g_a2 = Wide::Zero();
                if (row + 1 < size) {
                    g_a2 = -lambda[row] * solution[row + 1];
                }
                w.d_homogeneous[az].xp[layer] +=
                    g_e * w.homogeneous[az].omega[layer];
                w.d_homogeneous[az].omega[layer] +=
                    g_e * w.homogeneous[az].xp[layer];
                w.d_homogeneous[az].xm[layer] += g_c2;
                w.d_homogeneous[az].xp[layer + 1] -= g_d2;
                w.d_homogeneous[az].xm[layer + 1] -=
                    g_a2 * w.homogeneous[az].omega[layer + 1];
                w.d_homogeneous[az].omega[layer + 1] -=
                    g_a2 * w.homogeneous[az].xm[layer + 1];
            }

            const int layer = n - 1;
            const Wide g_c = -lambda[last_row] * solution[last_row - 1];
            const Wide g_d = -lambda[last_row] * solution[last_row];
            const Wide c_base =
                w.homogeneous[az].xm[layer] -
                2.0 * mu * albedo * delta * w.homogeneous[az].xp[layer];
            w.d_homogeneous[az].xm[layer] +=
                g_c * w.homogeneous[az].omega[layer];
            w.d_homogeneous[az].xp[layer] -= g_c * 2.0 * mu * albedo * delta *
                                             w.homogeneous[az].omega[layer];
            w.d_homogeneous[az].omega[layer] += g_c * c_base;
            d_albedo -= g_c * 2.0 * mu * delta * w.homogeneous[az].xp[layer] *
                        w.homogeneous[az].omega[layer];
            w.d_homogeneous[az].xp[layer] += g_d;
            w.d_homogeneous[az].xm[layer] -= g_d * 2.0 * mu * albedo * delta;
            d_albedo -= g_d * 2.0 * mu * delta * w.homogeneous[az].xm[layer];
        }
        return {d_albedo, d_thermal_surface};
    }

    template <bool Solar>
    void map_explicit_to_levels(const ColumnGeometry& geometry,
                                ExplicitWorkspace& w) {
        const int n = static_cast<int>(geometry.layer_thickness.size());
        if constexpr (Solar) {
            std::fill(w.attenuation.begin(), w.attenuation.end(), Wide::Zero());
            for (int boundary = 1; boundary <= n; ++boundary) {
                w.attenuation[boundary] =
                    w.d_transmission[boundary] * w.transmission[boundary];
            }
            for (int layer = 0; layer < n; ++layer) {
                const Wide inverse_od = positive_reciprocal(w.od[layer]);
                const Wide weight = w.d_secant[layer] * inverse_od;
                w.attenuation[layer] += weight;
                w.attenuation[layer + 1] -= weight;
                w.d_od[layer] -=
                    w.d_secant[layer] * w.secant[layer] * inverse_od;
            }
            for (int layer = 0; layer < n; ++layer) {
                Wide contribution = Wide::Zero();
                for (int boundary = 1; boundary <= n; ++boundary) {
                    contribution += w.attenuation[boundary] *
                                    geometry.chapman(boundary - 1, layer);
                }
                w.d_od[layer] -= contribution;
            }
        }

        for (int layer = 0; layer < n; ++layer) {
            if constexpr (!Solar) {
                Wide d_top = Wide::Zero();
                Wide d_bottom = Wide::Zero();
                Wide d_od = Wide::Zero();
                for (int lane = 0; lane < LANES; ++lane) {
                    const double top = w.level_emission[layer][lane];
                    const double bottom = w.level_emission[layer + 1][lane];
                    const double od = w.od[layer][lane];
                    if (od > THERMAL_MIN_OPTICAL_DEPTH &&
                        top > THERMAL_MIN_EMISSION &&
                        bottom > THERMAL_MIN_EMISSION) {
                        const double seed = w.d_thermal_b1[layer][lane];
                        d_top[lane] = seed / (top * od);
                        d_bottom[lane] = -seed / (bottom * od);
                        d_od[lane] = -seed * w.thermal_b1[layer][lane] / od;
                    }
                }
                w.d_level_emission[layer] += w.d_thermal_b0[layer] + d_top;
                w.d_level_emission[layer + 1] += d_bottom;
                w.d_od[layer] += d_od;
            }

            const double dz = geometry.layer_thickness[layer];
            const Wide od_density = w.od[layer] / dz;
            const Wide scattering_density = w.ssa[layer] * od_density;
            const Wide inverse_od_density = positive_reciprocal(od_density);
            const Wide inverse_scattering_density =
                positive_reciprocal(scattering_density);
            Wide layer_d_ssa = w.d_ssa[layer];
            for (int lane = 0; lane < LANES; ++lane) {
                if (w.ssa[layer][lane] >= 1.0 - 1.0e-9) {
                    layer_d_ssa[lane] = 0.0;
                }
            }
            const Wide common_od = 0.5 * w.d_od[layer] * dz;
            const Wide common_ssa = 0.5 * layer_d_ssa * inverse_od_density;
            const Wide common_b1 =
                0.5 * w.d_b1[layer] * inverse_scattering_density;
            for (int endpoint = 0; endpoint < 2; ++endpoint) {
                const int level = layer + endpoint;
                const Wide delta_b1 = w.level_b1[level] - w.b1[layer];
                w.d_level_extinction[level] +=
                    common_od +
                    common_ssa * (w.level_ssa[level] - w.ssa[layer]) +
                    common_b1 * w.level_ssa[level] * delta_b1;
                w.d_level_ssa[level] +=
                    common_ssa * w.level_extinction[level] +
                    common_b1 * w.level_extinction[level] * delta_b1;
                w.d_level_b1[level] +=
                    common_b1 * w.level_ssa[level] * w.level_extinction[level];
            }
        }
    }

    template <bool Solar>
    void reverse_explicit_layers(const ColumnGeometry& geometry,
                                 ExplicitWorkspace& w) {
        const int n = static_cast<int>(geometry.layer_thickness.size());
        constexpr int naz = Solar ? 2 : 1;
        const double mu = geometry.quadrature_cosine;
        const double solar_sine = std::sqrt(
            std::max(0.0, (1.0 - mu * mu) * (1.0 - geometry.solar_cosine *
                                                       geometry.solar_cosine)));
        for (int az = 0; az < naz; ++az) {
            for (int layer = 0; layer < n; ++layer) {
                const auto& h = w.homogeneous[az];
                const auto& p = w.particular[az];
                Wide d_xp = w.d_homogeneous[az].xp[layer];
                Wide d_xm = w.d_homogeneous[az].xm[layer];
                Wide d_omega = w.d_homogeneous[az].omega[layer];
                Wide d_k = w.d_homogeneous[az].k[layer];
                Wide d_norm = Wide::Zero();
                Wide d_cp = Wide::Zero();
                Wide d_cm = Wide::Zero();

                if constexpr (Solar) {
                    Wide d_ap = w.d_particular[az].ap[layer];
                    Wide d_am = w.d_particular[az].am[layer];
                    const Wide& d_gpt = w.d_particular[az].gpt[layer];
                    const Wide& d_gpb = w.d_particular[az].gpb[layer];
                    const Wide& d_gmt = w.d_particular[az].gmt[layer];
                    const Wide& d_gmb = w.d_particular[az].gmb[layer];
                    d_am += d_gpt * p.cm[layer] * h.xm[layer] +
                            d_gmt * p.cm[layer] * h.xp[layer];
                    d_cm += d_gpt * p.am[layer] * h.xm[layer] +
                            d_gmt * p.am[layer] * h.xp[layer];
                    d_xm += d_gpt * p.am[layer] * p.cm[layer];
                    d_xp += d_gmt * p.am[layer] * p.cm[layer];
                    d_ap += d_gpb * p.cp[layer] * h.xp[layer] +
                            d_gmb * p.cp[layer] * h.xm[layer];
                    d_cp += d_gpb * p.ap[layer] * h.xp[layer] +
                            d_gmb * p.ap[layer] * h.xm[layer];
                    d_xp += d_gpb * p.ap[layer] * p.cp[layer];
                    d_xm += d_gmb * p.ap[layer] * p.cp[layer];

                    Wide d_qp = Wide::Zero();
                    Wide d_qm = Wide::Zero();
                    const Wide inverse_norm = 1.0 / h.norm[layer];
                    const Wide d_ap_numerator = d_ap * inverse_norm;
                    d_norm -= d_ap * p.ap[layer] * inverse_norm;
                    d_qp += d_ap_numerator * h.xp[layer];
                    d_xp += d_ap_numerator * p.qp[layer];
                    d_qm += d_ap_numerator * h.xm[layer];
                    d_xm += d_ap_numerator * p.qm[layer];
                    const Wide d_am_numerator = d_am * inverse_norm;
                    d_norm -= d_am * p.am[layer] * inverse_norm;
                    d_qm += d_am_numerator * h.xp[layer];
                    d_xp += d_am_numerator * p.qm[layer];
                    d_qp += d_am_numerator * h.xm[layer];
                    d_xm += d_am_numerator * p.qp[layer];

                    if (az == 0) {
                        w.d_ssa[layer] +=
                            d_qp *
                                (1.0 +
                                 w.b1[layer] * geometry.solar_cosine * mu) /
                                FOUR_PI +
                            d_qm *
                                (1.0 -
                                 w.b1[layer] * geometry.solar_cosine * mu) /
                                FOUR_PI;
                        w.d_b1[layer] +=
                            d_qp * w.ssa[layer] * geometry.solar_cosine * mu /
                                FOUR_PI -
                            d_qm * w.ssa[layer] * geometry.solar_cosine * mu /
                                FOUR_PI;
                    } else {
                        w.d_ssa[layer] +=
                            (d_qp + d_qm) * w.b1[layer] * solar_sine / FOUR_PI;
                        w.d_b1[layer] +=
                            (d_qp + d_qm) * w.ssa[layer] * solar_sine / FOUR_PI;
                    }
                } else {
                    Wide d_at = w.d_particular[az].at[layer];
                    const Wide& d_gpt = w.d_particular[az].gpt[layer];
                    const Wide& d_gpb = w.d_particular[az].gpb[layer];
                    const Wide& d_gmt = w.d_particular[az].gmt[layer];
                    const Wide& d_gmb = w.d_particular[az].gmb[layer];
                    d_at += d_gpt * p.cm[layer] * h.xm[layer] +
                            d_gpb * p.cp[layer] * h.xp[layer] +
                            d_gmt * p.cm[layer] * h.xp[layer] +
                            d_gmb * p.cp[layer] * h.xm[layer];
                    d_cm += d_gpt * p.at[layer] * h.xm[layer] +
                            d_gmt * p.at[layer] * h.xp[layer];
                    d_cp += d_gpb * p.at[layer] * h.xp[layer] +
                            d_gmb * p.at[layer] * h.xm[layer];
                    d_xm += d_gpt * p.at[layer] * p.cm[layer] +
                            d_gmb * p.at[layer] * p.cp[layer];
                    d_xp += d_gpb * p.at[layer] * p.cp[layer] +
                            d_gmt * p.at[layer] * p.cm[layer];
                    const Wide inverse_norm = 1.0 / h.norm[layer];
                    const Wide d_at_numerator = d_at * inverse_norm;
                    d_norm -= d_at * p.at[layer] * inverse_norm;
                    w.d_ssa[layer] -=
                        d_at_numerator * (h.xp[layer] + h.xm[layer]);
                    d_xp += d_at_numerator * (1.0 - w.ssa[layer]);
                    d_xm += d_at_numerator * (1.0 - w.ssa[layer]);
                }

                const Wide& rate =
                    Solar ? w.secant[layer] : w.thermal_b1[layer];
                const Wide& amplitude =
                    Solar ? w.transmission[layer] : w.thermal_b0[layer];
                const Wide cp_denominator = rate - h.k[layer];
                const Wide inverse_cp = 1.0 / cp_denominator;
                const Wide cp_numerator = h.omega[layer] - p.exponential[layer];
                if constexpr (Solar) {
                    w.d_transmission[layer] += d_cp * cp_numerator * inverse_cp;
                } else {
                    w.d_thermal_b0[layer] += d_cp * cp_numerator * inverse_cp;
                }
                const Wide d_cp_numerator = d_cp * amplitude * inverse_cp;
                const Wide d_cp_denominator = -d_cp * p.cp[layer] * inverse_cp;
                d_omega += d_cp_numerator;
                Wide d_exponential = -d_cp_numerator;
                if constexpr (Solar) {
                    w.d_secant[layer] += d_cp_denominator;
                } else {
                    w.d_thermal_b1[layer] += d_cp_denominator;
                }
                d_k -= d_cp_denominator;

                const Wide cm_denominator = rate + h.k[layer];
                const Wide inverse_cm = 1.0 / cm_denominator;
                const Wide cm_numerator =
                    1.0 - h.omega[layer] * p.exponential[layer];
                if constexpr (Solar) {
                    w.d_transmission[layer] += d_cm * cm_numerator * inverse_cm;
                } else {
                    w.d_thermal_b0[layer] += d_cm * cm_numerator * inverse_cm;
                }
                const Wide d_cm_numerator = d_cm * amplitude * inverse_cm;
                const Wide d_cm_denominator = -d_cm * p.cm[layer] * inverse_cm;
                d_omega -= d_cm_numerator * p.exponential[layer];
                d_exponential -= d_cm_numerator * h.omega[layer];
                if constexpr (Solar) {
                    w.d_secant[layer] += d_cm_denominator;
                    w.d_secant[layer] -=
                        d_exponential * p.exponential[layer] * w.od[layer];
                } else {
                    w.d_thermal_b1[layer] += d_cm_denominator;
                    w.d_thermal_b1[layer] -=
                        d_exponential * p.exponential[layer] * w.od[layer];
                }
                d_k += d_cm_denominator;
                w.d_od[layer] -= d_exponential * p.exponential[layer] * rate;

                d_xp += d_norm * 2.0 * mu * h.xp[layer];
                d_xm -= d_norm * 2.0 * mu * h.xm[layer];
                w.d_od[layer] -= d_omega * h.omega[layer] * h.k[layer];
                d_k -= d_omega * h.omega[layer] * w.od[layer];
                const Wide inverse_k = 1.0 / h.k[layer];
                Wide d_s = (-d_xp + d_xm) * 0.5 * inverse_k;
                d_k += (d_xp - d_xm) * 0.5 * h.s[layer] * inverse_k.square();
                d_s += d_k * 0.5 * h.d[layer] * inverse_k;
                const Wide d_d = d_k * 0.5 * h.s[layer] * inverse_k;
                if (az == 0) {
                    w.d_ssa[layer] += d_d * w.b1[layer] * mu + d_s / mu;
                    w.d_b1[layer] += d_d * w.ssa[layer] * mu;
                } else {
                    w.d_ssa[layer] +=
                        d_s * w.b1[layer] * (1.0 - mu * mu) / (2.0 * mu);
                    w.d_b1[layer] +=
                        d_s * w.ssa[layer] * (1.0 - mu * mu) / (2.0 * mu);
                }
            }
        }
    }

    void resize_solution(ColumnSolution& solution, int n, int naz,
                         const Var& zero) {
        resize_vars(solution.od, n, zero);
        resize_vars(solution.ssa, n, zero);
        resize_vars(solution.b1, n, zero);
        resize_vars(solution.transmission, n + 1, zero);
        resize_vars(solution.secant, n, zero);
        resize_vars(solution.thermal_b0, n, zero);
        resize_vars(solution.thermal_b1, n, zero);
        for (int az = 0; az < naz; ++az) {
            auto& h = solution.homogeneous[az];
            resize_vars(h.k, n, zero);
            resize_vars(h.xp, n, zero);
            resize_vars(h.xm, n, zero);
            resize_vars(h.omega, n, zero);
            auto& p = solution.particular[az];
            resize_vars(p.ap, n, zero);
            resize_vars(p.am, n, zero);
            resize_vars(p.at, n, zero);
            resize_vars(p.gpt, n, zero);
            resize_vars(p.gpb, n, zero);
            resize_vars(p.gmt, n, zero);
            resize_vars(p.gmb, n, zero);
            auto& q = solution.bvp[az];
            const int size = 2 * n;
            resize_vars(q.e, size, zero);
            resize_vars(q.c, size, zero);
            resize_vars(q.d, size, zero);
            resize_vars(q.a, size, zero);
            resize_vars(q.b, size, zero);
            resize_vars(q.rhs, size, zero);
            resize_vars(q.alpha, size, zero);
            resize_vars(q.beta, size, zero);
            resize_vars(q.gamma, size, zero);
            resize_vars(q.inverse_mu, size, zero);
            resize_vars(q.z, size, zero);
        }
    }

    void pentadiagonal_solve(Bvp& q) {
        const int n = static_cast<int>(q.d.size());
        q.inverse_mu[0] = 1.0 / q.d[0];
        q.alpha[0] = q.a[0] * q.inverse_mu[0];
        q.beta[0] = q.b[0] * q.inverse_mu[0];
        q.z[0] = q.rhs[0] * q.inverse_mu[0];
        if (n > 1) {
            q.gamma[1] = q.c[1];
            q.inverse_mu[1] = 1.0 / (q.d[1] - q.alpha[0] * q.gamma[1]);
            q.alpha[1] = (q.a[1] - q.beta[0] * q.gamma[1]) * q.inverse_mu[1];
            q.beta[1] = q.b[1] * q.inverse_mu[1];
            q.z[1] = (q.rhs[1] - q.z[0] * q.gamma[1]) * q.inverse_mu[1];
        }
        for (int index = 2; index < n; ++index) {
            q.gamma[index] = q.c[index] - q.alpha[index - 2] * q.e[index];
            q.inverse_mu[index] =
                1.0 / (q.d[index] - q.beta[index - 2] * q.e[index] -
                       q.alpha[index - 1] * q.gamma[index]);
            if (index + 1 < n) {
                q.alpha[index] =
                    (q.a[index] - q.beta[index - 1] * q.gamma[index]) *
                    q.inverse_mu[index];
            }
            if (index + 2 < n) {
                q.beta[index] = q.b[index] * q.inverse_mu[index];
            }
            q.z[index] = (q.rhs[index] - q.z[index - 2] * q.e[index] -
                          q.z[index - 1] * q.gamma[index]) *
                         q.inverse_mu[index];
        }
        q.rhs[n - 1] = q.z[n - 1];
        if (n > 1) {
            q.rhs[n - 2] = q.z[n - 2] - q.alpha[n - 2] * q.rhs[n - 1];
        }
        for (int index = n - 3; index >= 0; --index) {
            q.rhs[index] = q.z[index] - q.alpha[index] * q.rhs[index + 1] -
                           q.beta[index] * q.rhs[index + 2];
        }
    }

    template <bool Solar>
    void prepare_and_solve_column(const ColumnGeometry& geometry,
                                  const std::vector<Var>& extinction,
                                  const std::vector<Var>& level_ssa,
                                  const std::vector<Var>& level_b1,
                                  const std::vector<Var>& emission,
                                  const Var& irradiance, const Var& albedo,
                                  const Var& thermal_surface,
                                  ColumnSolution& solution, const Var& zero) {
        const int n = static_cast<int>(geometry.layer_thickness.size());
        constexpr int naz = Solar ? 2 : 1;
        resize_solution(solution, n, naz, zero);

        for (int layer = 0; layer < n; ++layer) {
            const Var scattering_top = extinction[layer] * level_ssa[layer];
            const Var scattering_bottom =
                extinction[layer + 1] * level_ssa[layer + 1];
            const Var average_extinction =
                0.5 * (extinction[layer] + extinction[layer + 1]);
            const Var average_scattering =
                0.5 * (scattering_top + scattering_bottom);
            solution.od[layer] =
                average_extinction * geometry.layer_thickness[layer];
            solution.ssa[layer] = clamp_ssa(
                positive_ratio(average_scattering, average_extinction));
            solution.b1[layer] =
                positive_ratio(0.5 * (scattering_top * level_b1[layer] +
                                      scattering_bottom * level_b1[layer + 1]),
                               average_scattering);
            if constexpr (!Solar) {
                solution.thermal_b0[layer] = emission[layer];
                solution.thermal_b1[layer] = thermal_profile_slope(
                    emission[layer], emission[layer + 1], solution.od[layer]);
            }
        }

        if constexpr (Solar) {
            solution.transmission[0] = irradiance;
            std::vector<Var> slant(n, zero);
            for (int boundary = 0; boundary < n; ++boundary) {
                Var optical_depth = zero;
                for (int layer = 0; layer < n; ++layer) {
                    const double factor = geometry.chapman(boundary, layer);
                    if (factor != 0.0) {
                        optical_depth =
                            optical_depth + solution.od[layer] * factor;
                    }
                }
                slant[boundary] = optical_depth;
                solution.transmission[boundary + 1] =
                    (-optical_depth).exp() * irradiance;
            }
            for (int layer = 0; layer < n; ++layer) {
                const Var top_slant = layer == 0 ? zero : slant[layer - 1];
                solution.secant[layer] = positive_ratio(
                    slant[layer] - top_slant, solution.od[layer]);
            }
        }

        const double mu = geometry.quadrature_cosine;
        for (int az = 0; az < naz; ++az) {
            auto& h = solution.homogeneous[az];
            auto& p = solution.particular[az];
            for (int layer = 0; layer < n; ++layer) {
                Var d = zero;
                Var s = zero;
                if (az == 0) {
                    d = solution.ssa[layer] * solution.b1[layer] * mu -
                        1.0 / mu;
                    s = (solution.ssa[layer] - 1.0) / mu;
                } else {
                    d = zero - 1.0 / mu;
                    s = (solution.ssa[layer] * solution.b1[layer] *
                             (1.0 - mu * mu) -
                         2.0) /
                        (2.0 * mu);
                }
                const Var k = (s * d).sqrt();
                const Var xp = 0.5 * (1.0 - s / k);
                const Var xm = 0.5 * (1.0 + s / k);
                const Var omega = (-k * solution.od[layer]).exp();
                h.k[layer] = k;
                h.xp[layer] = xp;
                h.xm[layer] = xm;
                h.omega[layer] = omega;
                const Var norm = mu * (xp * xp - xm * xm);

                if constexpr (Solar) {
                    Var qp = zero;
                    Var qm = zero;
                    if (az == 0) {
                        qp = solution.ssa[layer] *
                             (1.0 +
                              solution.b1[layer] * geometry.solar_cosine * mu) /
                             FOUR_PI;
                        qm = solution.ssa[layer] *
                             (1.0 -
                              solution.b1[layer] * geometry.solar_cosine * mu) /
                             FOUR_PI;
                    } else {
                        const double angular = std::sqrt(std::max(
                            0.0, (1.0 - mu * mu) *
                                     (1.0 - geometry.solar_cosine *
                                                geometry.solar_cosine)));
                        qp = solution.ssa[layer] * solution.b1[layer] *
                             angular / FOUR_PI;
                        qm = qp;
                    }
                    const Var ap = (qp * xp + qm * xm) / norm;
                    const Var am = (qm * xp + qp * xm) / norm;
                    const Var cp =
                        solution.transmission[layer] *
                        stable_exp_difference(k, solution.secant[layer],
                                              solution.od[layer]);
                    const Var cm =
                        solution.transmission[layer] *
                        stable_exp_difference(zero, solution.secant[layer] + k,
                                              solution.od[layer]);
                    p.ap[layer] = ap;
                    p.am[layer] = am;
                    p.gpt[layer] = am * cm * xm;
                    p.gpb[layer] = ap * cp * xp;
                    p.gmt[layer] = am * cm * xp;
                    p.gmb[layer] = ap * cp * xm;
                } else {
                    const Var at =
                        (1.0 - solution.ssa[layer]) * (xp + xm) / norm;
                    const Var cp =
                        solution.thermal_b0[layer] *
                        stable_exp_difference(k, solution.thermal_b1[layer],
                                              solution.od[layer]);
                    const Var cm = solution.thermal_b0[layer] *
                                   stable_exp_difference(
                                       zero, solution.thermal_b1[layer] + k,
                                       solution.od[layer]);
                    p.at[layer] = at;
                    p.gpt[layer] = at * cm * xm;
                    p.gpb[layer] = at * cp * xp;
                    p.gmt[layer] = at * cm * xp;
                    p.gmb[layer] = at * cp * xm;
                }
            }
        }

        for (int az = 0; az < naz; ++az) {
            const auto& h = solution.homogeneous[az];
            const auto& p = solution.particular[az];
            auto& q = solution.bvp[az];
            q.rhs[0] = -p.gpt[0];
            for (int layer = 0; layer < n - 1; ++layer) {
                q.rhs[2 * layer + 1] = p.gmt[layer + 1] - p.gmb[layer];
                q.rhs[2 * layer + 2] = p.gpt[layer + 1] - p.gpb[layer];
            }
            const int last = 2 * n - 1;
            const double delta = az == 0 ? 1.0 : 0.0;
            Var direct_boundary_source = zero;
            if constexpr (Solar) {
                direct_boundary_source = delta * geometry.solar_cosine *
                                         albedo / EIGEN_PI *
                                         solution.transmission[n];
            } else {
                direct_boundary_source = thermal_surface;
            }
            q.rhs[last] =
                direct_boundary_source -
                (p.gmb[n - 1] - 2.0 * delta * mu * albedo * p.gpb[n - 1]);
            q.d[0] = h.xp[0];
            q.a[0] = h.xm[0] * h.omega[0];
            for (int layer = 0; layer < n - 1; ++layer) {
                const int row = 2 * layer;
                q.c[row + 1] = h.xm[layer] * h.omega[layer];
                q.d[row + 1] = h.xp[layer];
                q.a[row + 1] = -h.xm[layer + 1];
                q.b[row + 1] = -h.xp[layer + 1] * h.omega[layer + 1];
                q.e[row + 2] = h.xp[layer] * h.omega[layer];
                q.c[row + 2] = h.xm[layer];
                q.d[row + 2] = -h.xp[layer + 1];
                q.a[row + 2] = -h.xm[layer + 1] * h.omega[layer + 1];
            }
            q.c[last] =
                (h.xm[n - 1] - 2.0 * mu * albedo * delta * h.xp[n - 1]) *
                h.omega[n - 1];
            q.d[last] = h.xp[n - 1] - 2.0 * mu * albedo * delta * h.xm[n - 1];
            pentadiagonal_solve(q);
        }
    }
} // namespace
template <sasktran2::twostream::SourceType SOURCE_TYPE>
void CppTwoStreamSourceAdapter<SOURCE_TYPE>::initialize_geometry(
    const sasktran2::viewinggeometry::InternalViewingGeometry&
        internal_viewing) {
    ZoneScopedN("C++ SIMD Twostream Initialize Geometry");
    auto& impl = *m_impl;
    if (impl.config == nullptr) {
        throw std::logic_error(
            "C++ two-stream source configuration was not initialized");
    }
    const int nlyr = impl.geometry.size() - 1;
    const bool is_spherical = impl.geometry.coordinates().geometry_type() ==
                              sasktran2::geometrytype::spherical;
    const double reference_cos_sza =
        impl.geometry.coordinates().cos_sza_at_reference();
    std::vector<double> sza_grid;
    if constexpr (sasktran2::twostream::has_solar<SOURCE_TYPE>()) {
        const int requested_sza = impl.config->num_do_sza();
        if (is_spherical && requested_sza > 1) {
            const auto bounds =
                sasktran2::raytracing::min_max_cos_sza_of_all_rays(
                    internal_viewing.traced_rays);
            if (bounds.first <= bounds.second &&
                bounds.second - bounds.first > 1.0e-14) {
                sza_grid.resize(requested_sza);
                for (int index = 0; index < requested_sza; ++index) {
                    sza_grid[index] =
                        bounds.first +
                        (bounds.second - bounds.first) * index /
                            static_cast<double>(requested_sza - 1);
                }
            }
        }
    }
    if (sza_grid.empty()) {
        sza_grid.push_back(reference_cos_sza);
    }

    impl.pconfigs.clear();
    impl.geometry_layers.clear();
    impl.columns.clear();
    impl.pconfigs.reserve(sza_grid.size());
    impl.geometry_layers.reserve(sza_grid.size());
    impl.columns.reserve(sza_grid.size());
    for (double cos_sza : sza_grid) {
        auto pconfig =
            std::make_unique<sasktran_disco::PersistentConfiguration<1>>();
        pconfig->configure(impl.spec, *impl.config, cos_sza, nlyr,
                           internal_viewing.traced_rays);
        auto layers = std::make_unique<sasktran_disco::GeometryLayerArray<1>>(
            *pconfig, impl.geometry);
        ColumnGeometry column;
        column.layer_thickness.resize(nlyr);
        column.chapman = layers->chapman_factors();
        column.solar_cosine = cos_sza;
        for (int layer = 0; layer < nlyr; ++layer) {
            column.layer_thickness[layer] =
                layers->layer_ceiling()(layer) - layers->layer_floor()(layer);
        }
        impl.pconfigs.push_back(std::move(pconfig));
        impl.geometry_layers.push_back(std::move(layers));
        impl.columns.push_back(std::move(column));
    }

    impl.views.clear();
    impl.spherical.reset();
    if (!is_spherical) {
        impl.views.reserve(internal_viewing.traced_rays.size());
        for (const auto& ray : internal_viewing.traced_rays) {
            const double viewing_cosine = -ray.observer_and_look.look_away.z();
            if (viewing_cosine <= 0.0) {
                throw sasktran_disco::InternalRuntimeError(
                    "C++ two-stream currently only supports upwelling "
                    "plane-parallel radiances");
            }
            impl.views.push_back(
                {viewing_cosine, -ray.observer_and_look.relative_azimuth});
        }
    } else {
        impl.spherical = std::make_unique<SphericalGeometry>();
        auto& spherical = *impl.spherical;
        spherical.sza_grid = sza_grid;
        spherical.columns = impl.columns;
        const auto& altitude_grid = impl.geometry.altitude_grid().grid();
        const std::size_t nlevel = altitude_grid.size();
        const double earth_radius = impl.geometry.coordinates().earth_radius();
        spherical.ray_offsets.reserve(internal_viewing.traced_rays.size() + 1);
        spherical.ground_hit.reserve(internal_viewing.traced_rays.size());
        spherical.ground_cos_sza.reserve(internal_viewing.traced_rays.size());
        spherical.ray_offsets.push_back(0);
        spherical.od_offsets.push_back(0);

        for (const auto& ray : internal_viewing.traced_rays) {
            spherical.ground_hit.push_back(ray.ground_is_hit ? 1 : 0);
            spherical.ground_cos_sza.push_back(
                ray.layers.empty() ? reference_cos_sza
                                   : ray.layers.front().cos_sza_exit);
            for (std::size_t layer_index = 0; layer_index < ray.layers.size();
                 ++layer_index) {
                const auto& layer = ray.layers[layer_index];
                const auto source_weights = normalized_source_weights(layer);
                const Eigen::Vector3d position =
                    source_weights.first * layer.entrance.position +
                    source_weights.second * layer.exit.position;
                const double altitude = position.norm() - earth_radius;
                const double* upper =
                    std::upper_bound(altitude_grid.data(),
                                     altitude_grid.data() + nlevel, altitude);
                std::size_t bottom_layer =
                    upper == altitude_grid.data()
                        ? 0
                        : static_cast<std::size_t>(upper -
                                                   altitude_grid.data() - 1);
                bottom_layer = std::min(bottom_layer, nlevel - 2);
                const double bottom = altitude_grid(bottom_layer);
                const double top = altitude_grid(bottom_layer + 1);
                const double fraction_from_top =
                    std::clamp((top - altitude) / (top - bottom), 0.0, 1.0);
                const Eigen::Vector3d local_up = position.normalized();
                const double cosine = std::clamp(
                    -layer.average_look_away.dot(local_up), -1.0, 1.0);
                const double sin_azimuth =
                    source_weights.first * std::sin(layer.saz_entrance) +
                    source_weights.second * std::sin(layer.saz_exit);
                const double cos_azimuth =
                    source_weights.first * std::cos(layer.saz_entrance) +
                    source_weights.second * std::cos(layer.saz_exit);

                spherical.segment_layers.push_back(
                    static_cast<std::size_t>(nlyr - 1) - bottom_layer);
                spherical.segment_fractions.push_back(fraction_from_top);
                spherical.segment_cosines.push_back(cosine);
                spherical.segment_relative_azimuths.push_back(
                    std::atan2(sin_azimuth, cos_azimuth));
                spherical.segment_cos_sza.push_back(
                    std::clamp(source_weights.first * layer.cos_sza_entrance +
                                   source_weights.second * layer.cos_sza_exit,
                               -1.0, 1.0));

                const auto stencil = ray.optical_depth_weights(layer_index);
                for (std::size_t entry = 0; entry < stencil.size(); ++entry) {
                    const auto [index, weight] = stencil[entry];
                    if (index < 0 ||
                        static_cast<std::size_t>(index) >= nlevel) {
                        throw sasktran_disco::InternalRuntimeError(
                            "Traced ray optical-depth stencil is outside "
                            "the one-dimensional atmosphere grid");
                    }
                    if (weight != 0.0) {
                        spherical.od_indices.push_back(nlevel - 1 - index);
                        spherical.od_weights.push_back(weight);
                    }
                }
                spherical.od_offsets.push_back(spherical.od_indices.size());
            }
            spherical.ray_offsets.push_back(spherical.segment_layers.size());
        }
    }

    impl.workers.resize(std::max(1, impl.config->num_wavelength_threads()));
    for (auto& worker : impl.workers) {
        impl.resize_worker(worker);
    }
}

template <sasktran2::twostream::SourceType SOURCE_TYPE>
void CppTwoStreamSourceAdapter<SOURCE_TYPE>::initialize_atmosphere(
    const sasktran2::atmosphere::Atmosphere<1>& atmosphere) {
    m_impl->atmosphere = &atmosphere;
}

template <sasktran2::twostream::SourceType SOURCE_TYPE>
void CppTwoStreamSourceAdapter<SOURCE_TYPE>::set_wavelength_block_capacity(
    int block_capacity) {
    if (block_capacity < 1) {
        throw std::invalid_argument(
            "C++ two-stream wavelength block capacity must be positive");
    }
    m_impl->block_capacity = block_capacity;
    for (auto& worker : m_impl->workers) {
        m_impl->resize_worker(worker);
    }
}

template <sasktran2::twostream::SourceType SOURCE_TYPE>
void CppTwoStreamSourceAdapter<SOURCE_TYPE>::calculate(
    const sasktran2::WavelengthBlock<>& block, int threadidx) {
    ZoneScopedN("C++ SIMD Twostream Calculate Block");
    auto& impl = *m_impl;
    if (impl.atmosphere == nullptr || threadidx < 0 ||
        static_cast<std::size_t>(threadidx) >= impl.workers.size() ||
        block.count < 1 || block.count > impl.block_capacity ||
        block.start < 0 ||
        block.start + block.count > impl.atmosphere->num_wavel()) {
        throw std::invalid_argument("Invalid C++ two-stream wavelength block");
    }
    auto& worker = impl.workers[threadidx];
    for (int packet_base = 0; packet_base < block.count; packet_base += LANES) {
        impl.calculate_packet(block, packet_base, worker);
    }
}

template <sasktran2::twostream::SourceType SOURCE_TYPE>
void CppTwoStreamSourceAdapter<SOURCE_TYPE>::start_of_ray_source(
    const sasktran2::WavelengthBlock<>& block, int losidx, int wavel_threadidx,
    int, sasktran2::WavelengthBlockDual<1>& source) const {
    ZoneScopedN("C++ SIMD Twostream Cached Source");
    const auto& impl = *m_impl;
    if (wavel_threadidx < 0 ||
        static_cast<std::size_t>(wavel_threadidx) >= impl.workers.size() ||
        losidx < 0) {
        throw std::invalid_argument("Invalid C++ two-stream result index");
    }
    const auto& worker = impl.workers[wavel_threadidx];
    const std::size_t surface_offset =
        static_cast<std::size_t>(losidx) * impl.block_capacity;
    source.value.row(0).head(block.count) +=
        Eigen::Map<const Eigen::RowVectorXd>(
            worker.radiance.data() + surface_offset, block.count);
    if (impl.atmosphere->num_deriv() == 0) {
        return;
    }

    const int nlevel = impl.geometry.size();
    for (int rust_level = 0; rust_level < nlevel; ++rust_level) {
        const int cpp_level = nlevel - 1 - rust_level;
        const std::size_t offset =
            (static_cast<std::size_t>(losidx) * nlevel + rust_level) *
            impl.block_capacity;
        source.deriv.row(cpp_level).head(block.count) +=
            Eigen::Map<const Eigen::RowVectorXd>(
                worker.extinction.data() + offset, block.count);
        source.deriv.row(impl.atmosphere->ssa_deriv_start_index() + cpp_level)
            .head(block.count) += Eigen::Map<const Eigen::RowVectorXd>(
            worker.ssa.data() + offset, block.count);
        for (int group = 0;
             group < impl.atmosphere->num_scattering_deriv_groups(); ++group) {
            const int derivative_index =
                impl.atmosphere->scat_deriv_start_index() + group * nlevel +
                cpp_level;
            for (int lane = 0; lane < block.count; ++lane) {
                const int wavelength = block.wavelength(lane);
                double phase_derivative =
                    impl.atmosphere->storage().d_leg_coeff(1, cpp_level,
                                                           wavelength, group);
                if (impl.atmosphere->storage().applied_f_order > 0) {
                    const double delta_m =
                        impl.atmosphere->storage().f(cpp_level, wavelength);
                    phase_derivative -= 3.0 *
                                        impl.atmosphere->storage().d_f(
                                            cpp_level, wavelength, group) /
                                        ((1.0 - delta_m) * (1.0 - delta_m));
                }
                source.deriv(derivative_index, lane) +=
                    phase_derivative * worker.b1[offset + lane];
            }
        }
        if constexpr (sasktran2::twostream::has_thermal<SOURCE_TYPE>()) {
            if (impl.atmosphere->include_emission_derivatives()) {
                source.deriv
                    .row(impl.atmosphere->emission_deriv_start_index() +
                         cpp_level)
                    .head(block.count) += Eigen::Map<const Eigen::RowVectorXd>(
                    worker.emission.data() + offset, block.count);
            }
        }
    }

    for (int deriv = 0; deriv < impl.atmosphere->surface().num_deriv();
         ++deriv) {
        source.deriv.row(impl.atmosphere->surface_deriv_start_index() + deriv)
            .head(block.count) += Eigen::Map<const Eigen::RowVectorXd>(
            worker.surface_albedo.data() + surface_offset, block.count);
    }
    if constexpr (sasktran2::twostream::has_thermal<SOURCE_TYPE>()) {
        if (impl.atmosphere->include_emission_derivatives()) {
            source.deriv
                .row(impl.atmosphere->surface_emission_deriv_start_index())
                .head(block.count) += Eigen::Map<const Eigen::RowVectorXd>(
                worker.surface_emission.data() + surface_offset, block.count);
        }
    }
}

template class CppTwoStreamSourceAdapter<
    sasktran2::twostream::SourceType::ONLY_SOLAR>;
template class CppTwoStreamSourceAdapter<
    sasktran2::twostream::SourceType::ONLY_THERMAL>;
