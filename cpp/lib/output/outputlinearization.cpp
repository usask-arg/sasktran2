#include "sasktran2/output.h"

namespace sasktran2 {
    namespace {
        template <int NSTOKES, typename Radiance>
        void mapped_volume_derivatives(
            const Radiance& radiance,
            const sasktran2::DerivativeMapping& mapping, int wavelidx,
            int ngeometry, int num_scattering_groups,
            bool include_emission_derivatives, Eigen::MatrixXd& native) {
            const auto d_rad_by_ssa = radiance.d_ssa(ngeometry);
            const auto d_rad_by_extinction = radiance.d_extinction(ngeometry);
            const auto& d_ssa = mapping.native_mapping().d_ssa.value();
            const auto& d_extinction =
                mapping.native_mapping().d_extinction.value();

            for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                native.row(stokes).array() =
                    d_ssa.col(wavelidx).transpose().array() *
                        d_rad_by_ssa(stokes, Eigen::placeholders::all).array() +
                    d_extinction.col(wavelidx).transpose().array() *
                        d_rad_by_extinction(stokes, Eigen::placeholders::all)
                            .array();
            }

            if (mapping.is_scattering_derivative()) {
                const auto d_rad_by_scattering = radiance.d_scatterer(
                    ngeometry, mapping.get_scattering_index());
                const auto& scattering_factor =
                    mapping.native_mapping().scat_factor.value();
                for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                    native.row(stokes).array() +=
                        scattering_factor.col(wavelidx).transpose().array() *
                        d_rad_by_scattering(stokes, Eigen::placeholders::all)
                            .array();
                }
            }

            if (mapping.native_mapping().d_emission.has_value() &&
                include_emission_derivatives) {
                const auto d_rad_by_emission =
                    radiance.d_emission(ngeometry, num_scattering_groups);
                const auto& d_emission =
                    mapping.native_mapping().d_emission.value();
                native.row(0).array() +=
                    d_emission.col(wavelidx).transpose().array() *
                    d_rad_by_emission(0, Eigen::placeholders::all).array();
            }
        }

        template <int NSTOKES, typename Radiance>
        Eigen::Matrix<double, NSTOKES, 1> mapped_surface_derivative(
            const Radiance& radiance,
            const sasktran2::SurfaceDerivativeMapping& mapping, int wavelidx,
            int ngeometry, int num_source_groups, int num_surface_derivatives,
            int surface_emission_index, bool include_emission_derivatives) {
            Eigen::Matrix<double, NSTOKES, 1> result;
            result.setZero();
            if (mapping.native_surface_mapping().d_brdf.has_value()) {
                const auto d_rad_by_surface = radiance.d_brdf(
                    ngeometry, num_source_groups, num_surface_derivatives);
                for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                    result(stokes) =
                        d_rad_by_surface(stokes, Eigen::placeholders::all)
                            .dot(mapping.native_surface_mapping()
                                     .d_brdf.value()
                                     .row(wavelidx));
                }
            }
            if (mapping.native_surface_mapping().d_emission.has_value() &&
                include_emission_derivatives) {
                result(0) +=
                    radiance.deriv(0, surface_emission_index) *
                    mapping.native_surface_mapping().d_emission.value()(
                        wavelidx, 0);
            }
            return result;
        }

        template <int NSTOKES>
        Eigen::Matrix<double, NSTOKES, 1>
        rotate_stokes(const Eigen::Matrix<double, NSTOKES, 1>& value, double c,
                      double s) {
            auto result = value;
            if constexpr (NSTOKES >= 3) {
                result(1) = c * value(1) - s * value(2);
                result(2) = s * value(1) + c * value(2);
            }
            return result;
        }

        template <int NSTOKES>
        Eigen::Matrix<double, NSTOKES, 1>
        transpose_rotate_stokes(const Eigen::Matrix<double, NSTOKES, 1>& value,
                                double c, double s) {
            auto result = value;
            if constexpr (NSTOKES >= 3) {
                result(1) = c * value(1) + s * value(2);
                result(2) = -s * value(1) + c * value(2);
            }
            return result;
        }

        template <int NSTOKES, typename ParameterMap>
        void validate_volume_parameter_sizes(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
            const ParameterMap& parameters, const std::string& quantity) {
            const auto& mappings =
                atmosphere.storage().derivative_mappings_const();
            const int num_geometry =
                atmosphere.storage().total_extinction.rows();
            for (const auto& [name, parameter] : parameters) {
                const auto mapping = mappings.find(name);
                if (mapping == mappings.end()) {
                    throw std::invalid_argument(
                        "Unknown volume derivative mapping '" + name + "'");
                }

                int expected_size = num_geometry;
                if (mapping->second.get_interpolator_const().has_value()) {
                    const auto& interpolator =
                        mapping->second.get_interpolator_const().value();
                    if (interpolator.rows() != num_geometry) {
                        throw std::invalid_argument(
                            "Volume derivative mapping '" + name +
                            "' has an interpolator with " +
                            std::to_string(interpolator.rows()) +
                            " rows; expected " + std::to_string(num_geometry));
                    }
                    expected_size = interpolator.cols();
                }
                if (parameter.size() != expected_size) {
                    throw std::invalid_argument(
                        "Volume " + quantity + " for mapping '" + name +
                        "' has size " + std::to_string(parameter.size()) +
                        "; expected " + std::to_string(expected_size));
                }
            }
        }

        template <int NSTOKES, typename ParameterMap>
        void validate_surface_parameter_sizes(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
            const ParameterMap& parameters, const std::string& quantity) {
            const auto& mappings = atmosphere.surface().derivative_mappings();
            const int num_wavelengths = atmosphere.num_wavel();
            for (const auto& [name, parameter] : parameters) {
                const auto mapping = mappings.find(name);
                if (mapping == mappings.end()) {
                    throw std::invalid_argument(
                        "Unknown surface derivative mapping '" + name + "'");
                }

                int expected_size = num_wavelengths;
                if (mapping->second.get_interpolator_const().has_value()) {
                    const auto& interpolator =
                        mapping->second.get_interpolator_const().value();
                    if (interpolator.rows() != num_wavelengths) {
                        throw std::invalid_argument(
                            "Surface derivative mapping '" + name +
                            "' has an interpolator with " +
                            std::to_string(interpolator.rows()) +
                            " rows; expected " +
                            std::to_string(num_wavelengths));
                    }
                    expected_size = interpolator.cols();
                }
                if (parameter.size() != expected_size) {
                    throw std::invalid_argument(
                        "Surface " + quantity + " for mapping '" + name +
                        "' has size " + std::to_string(parameter.size()) +
                        "; expected " + std::to_string(expected_size));
                }
            }
        }

        template <typename Vector>
        void validate_product_output_buffer(const Vector& buffer,
                                            Eigen::Index expected_size,
                                            const std::string& quantity) {
            if (buffer.size() != expected_size ||
                (expected_size > 0 && buffer.data() == nullptr)) {
                throw std::invalid_argument(quantity + " buffer has size " +
                                            std::to_string(buffer.size()) +
                                            "; expected " +
                                            std::to_string(expected_size));
            }
        }
    } // namespace

    template <int NSTOKES> void OutputJVP<NSTOKES>::resize() {
        const Eigen::Index expected_size =
            static_cast<Eigen::Index>(NSTOKES) *
            static_cast<Eigen::Index>(this->m_nlos) *
            static_cast<Eigen::Index>(this->m_nwavel);
        validate_product_output_buffer(m_radiance, expected_size,
                                       "JVP radiance");
        validate_product_output_buffer(m_jvp, expected_size, "JVP product");
        if (this->m_atmosphere->num_deriv() == 0) {
            throw std::invalid_argument(
                "JVP calculation requires an atmosphere constructed with "
                "calculate_derivatives=true");
        }
        validate_volume_parameter_sizes(*this->m_atmosphere,
                                        m_derivative_tangents, "JVP tangent");
        validate_surface_parameter_sizes(*this->m_atmosphere,
                                         m_surface_tangents, "JVP tangent");

        m_jvp.setZero();
        m_native_parameter_tangents.clear();
        for (const auto& [name, parameter_tangent] : m_derivative_tangents) {
            const auto& mapping =
                this->m_atmosphere->storage().derivative_mappings_const().at(
                    name);
            if (mapping.get_interpolator_const().has_value()) {
                m_native_parameter_tangents[name] =
                    mapping.get_interpolator_const().value() *
                    parameter_tangent;
            } else {
                m_native_parameter_tangents[name] = parameter_tangent;
            }
        }
        m_native_thread_storage.resize(this->m_config->num_threads());
        m_mapped_thread_storage.resize(this->m_config->num_threads());
        int max_output = 1;
        for (const auto& [name, tangent] : m_derivative_tangents) {
            max_output = std::max(max_output, static_cast<int>(tangent.size()));
        }
        for (int thread = 0; thread < this->m_config->num_threads(); ++thread) {
            m_native_thread_storage[thread].resize(NSTOKES, this->m_ngeometry);
            m_mapped_thread_storage[thread].resize(NSTOKES, max_output);
        }
    }

    template <int NSTOKES>
    void OutputJVP<NSTOKES>::native_tangent(
        int wavelidx, Eigen::Ref<Eigen::VectorXd> tangent) const {
        if (tangent.size() != this->m_atmosphere->num_deriv()) {
            throw std::invalid_argument(
                "Native JVP tangent has incorrect size");
        }
        tangent.setZero();
        const int num_geometry = this->m_ngeometry;
        for (const auto& [name, unused] : m_derivative_tangents) {
            (void)unused;
            const auto& mapping =
                this->m_atmosphere->storage().derivative_mappings_const().at(
                    name);
            const auto& native_parameter = m_native_parameter_tangents.at(name);
            const auto& native_mapping = mapping.native_mapping();
            if (native_mapping.d_extinction.has_value()) {
                tangent.head(num_geometry).array() +=
                    native_mapping.d_extinction.value().col(wavelidx).array() *
                    native_parameter.array();
            }
            if (native_mapping.d_ssa.has_value()) {
                tangent
                    .segment(this->m_atmosphere->ssa_deriv_start_index(),
                             num_geometry)
                    .array() +=
                    native_mapping.d_ssa.value().col(wavelidx).array() *
                    native_parameter.array();
            }
            if (mapping.is_scattering_derivative()) {
                const int start = this->m_atmosphere->scat_deriv_start_index() +
                                  mapping.get_scattering_index() * num_geometry;
                tangent.segment(start, num_geometry).array() +=
                    native_mapping.scat_factor.value().col(wavelidx).array() *
                    native_parameter.array();
            }
            if (native_mapping.d_emission.has_value() &&
                this->m_atmosphere->include_emission_derivatives()) {
                tangent
                    .segment(this->m_atmosphere->emission_deriv_start_index(),
                             num_geometry)
                    .array() +=
                    native_mapping.d_emission.value().col(wavelidx).array() *
                    native_parameter.array();
            }
        }

        for (const auto& [name, parameter_tangent] : m_surface_tangents) {
            const auto& mapping =
                this->m_atmosphere->surface().derivative_mappings().at(name);
            double parameter_direction;
            if (mapping.get_interpolator_const().has_value()) {
                parameter_direction =
                    mapping.get_interpolator_const().value().row(wavelidx).dot(
                        parameter_tangent);
            } else {
                parameter_direction = parameter_tangent(wavelidx);
            }
            const auto& native_mapping = mapping.native_surface_mapping();
            if (native_mapping.d_brdf.has_value()) {
                tangent
                    .segment(this->m_atmosphere->surface_deriv_start_index(),
                             this->m_atmosphere->surface().num_deriv())
                    .array() +=
                    parameter_direction * native_mapping.d_brdf.value()
                                              .row(wavelidx)
                                              .transpose()
                                              .array();
            }
            if (native_mapping.d_emission.has_value() &&
                this->m_atmosphere->include_emission_derivatives()) {
                tangent(
                    this->m_atmosphere->surface_emission_deriv_start_index()) +=
                    parameter_direction *
                    native_mapping.d_emission.value()(wavelidx, 0);
            }
        }
    }

    template <int NSTOKES>
    void OutputJVP<NSTOKES>::assign_native(
        int losidx, int wavelidx, const Eigen::Vector<double, NSTOKES>& value,
        const Eigen::Vector<double, NSTOKES>& jvp) {
        const int linear_index =
            NSTOKES * this->m_nlos * wavelidx + NSTOKES * losidx;
        double c = 1.0;
        double s = 0.0;
        if constexpr (NSTOKES >= 3) {
            c = this->m_stokes_C[losidx];
            s = this->m_stokes_S[losidx];
        }
        const auto rotated_value = rotate_stokes<NSTOKES>(value, c, s);
        const auto rotated_jvp = rotate_stokes<NSTOKES>(jvp, c, s);
        for (int stokes = 0; stokes < NSTOKES; ++stokes) {
            m_radiance(linear_index + stokes) = rotated_value(stokes);
            m_jvp(linear_index + stokes) = rotated_jvp(stokes);
        }
    }

    template <int NSTOKES>
    template <typename Radiance>
    void OutputJVP<NSTOKES>::assign_lane(const Radiance& radiance, int losidx,
                                         int wavelidx, int threadidx) {
        const int linear_index =
            NSTOKES * this->m_nlos * wavelidx + NSTOKES * losidx;
        double c = 1.0;
        double s = 0.0;
        if constexpr (NSTOKES >= 3) {
            c = this->m_stokes_C[losidx];
            s = this->m_stokes_S[losidx];
        }

        Eigen::Matrix<double, NSTOKES, 1> value;
        Eigen::Matrix<double, NSTOKES, 1> directional;
        directional.setZero();
        for (int stokes = 0; stokes < NSTOKES; ++stokes) {
            value(stokes) = radiance.value(stokes);
        }

        auto& native = m_native_thread_storage[threadidx];
        auto& mapped = m_mapped_thread_storage[threadidx];
        for (const auto& [name, tangent] : m_derivative_tangents) {
            const auto& mapping =
                this->m_atmosphere->storage().derivative_mappings_const().at(
                    name);
            mapped_volume_derivatives<NSTOKES>(
                radiance, mapping, wavelidx, this->m_ngeometry,
                this->m_atmosphere->num_scattering_deriv_groups(),
                this->m_atmosphere->include_emission_derivatives(), native);
            if (mapping.get_interpolator_const().has_value()) {
                const auto& interpolator =
                    mapping.get_interpolator_const().value();
                if (tangent.size() != interpolator.cols()) {
                    throw std::invalid_argument(
                        "Volume linearization tangent has incorrect size");
                }
                mapped.leftCols(interpolator.cols()).noalias() =
                    native * interpolator;
                directional.noalias() +=
                    mapped.leftCols(interpolator.cols()) * tangent;
            } else {
                if (tangent.size() != native.cols()) {
                    throw std::invalid_argument(
                        "Volume linearization tangent has incorrect size");
                }
                directional.noalias() += native * tangent;
            }
        }

        for (const auto& [name, tangent] : m_surface_tangents) {
            const auto& mapping =
                this->m_atmosphere->surface().derivative_mappings().at(name);
            double parameter_direction = 0.0;
            if (mapping.get_interpolator_const().has_value()) {
                const auto& interpolator =
                    mapping.get_interpolator_const().value();
                if (tangent.size() != interpolator.cols()) {
                    throw std::invalid_argument(
                        "Surface linearization tangent has incorrect size");
                }
                parameter_direction = interpolator.row(wavelidx).dot(tangent);
            } else {
                if (tangent.size() != this->m_nwavel) {
                    throw std::invalid_argument(
                        "Native spectral surface tangent has incorrect size");
                }
                parameter_direction = tangent(wavelidx);
            }
            directional +=
                parameter_direction *
                mapped_surface_derivative<NSTOKES>(
                    radiance, mapping, wavelidx, this->m_ngeometry,
                    this->m_atmosphere->num_source_deriv_groups(),
                    this->m_atmosphere->surface().num_deriv(),
                    this->m_atmosphere->surface_emission_deriv_start_index(),
                    this->m_atmosphere->include_emission_derivatives());
        }

        value = rotate_stokes<NSTOKES>(value, c, s);
        directional = rotate_stokes<NSTOKES>(directional, c, s);
        for (int stokes = 0; stokes < NSTOKES; ++stokes) {
            m_radiance(linear_index + stokes) = value(stokes);
            m_jvp(linear_index + stokes) = directional(stokes);
        }
    }

    template <int NSTOKES>
    void OutputJVP<NSTOKES>::assign(
        const sasktran2::WavelengthBlock<>& block,
        const sasktran2::WavelengthBlockDual<NSTOKES>& radiance, int losidx,
        int threadidx) {
        for (int lane = 0; lane < block.count; ++lane) {
            const sasktran2::WavelengthBlockConstLaneDualView<NSTOKES>
                radiance_lane(radiance, lane);
            assign_lane(radiance_lane, losidx, block.start + lane, threadidx);
        }
    }

    template <int NSTOKES>
    void OutputVJP<NSTOKES>::set_derivative_gradient_memory(
        const std::string& name,
        Eigen::Map<Eigen::VectorXd> derivative_gradient) {
        m_derivative_gradients.insert(
            {name, Eigen::Map<Eigen::VectorXd>(nullptr, 0)});
        auto* target = &m_derivative_gradients.at(name);
        new (target) Eigen::Map<Eigen::VectorXd>(derivative_gradient.data(),
                                                 derivative_gradient.size());
    }

    template <int NSTOKES>
    void OutputVJP<NSTOKES>::set_surface_gradient_memory(
        const std::string& name,
        Eigen::Map<Eigen::VectorXd> derivative_gradient) {
        m_surface_gradients.insert(
            {name, Eigen::Map<Eigen::VectorXd>(nullptr, 0)});
        auto* target = &m_surface_gradients.at(name);
        new (target) Eigen::Map<Eigen::VectorXd>(derivative_gradient.data(),
                                                 derivative_gradient.size());
    }

    template <int NSTOKES> void OutputVJP<NSTOKES>::resize() {
        const Eigen::Index expected_size =
            static_cast<Eigen::Index>(NSTOKES) *
            static_cast<Eigen::Index>(this->m_nlos) *
            static_cast<Eigen::Index>(this->m_nwavel);
        validate_product_output_buffer(m_radiance, expected_size,
                                       "VJP radiance");
        validate_product_output_buffer(m_cotangent, expected_size,
                                       "VJP cotangent");
        if (this->m_atmosphere->num_deriv() == 0) {
            throw std::invalid_argument(
                "VJP calculation requires an atmosphere constructed with "
                "calculate_derivatives=true");
        }
        validate_volume_parameter_sizes(*this->m_atmosphere,
                                        m_derivative_gradients, "VJP gradient");
        validate_surface_parameter_sizes(*this->m_atmosphere,
                                         m_surface_gradients, "VJP gradient");

        const int num_threads = this->m_config->num_threads();
        m_native_thread_storage.resize(num_threads);
        m_thread_derivative_gradients.resize(num_threads);
        m_thread_surface_gradients.resize(num_threads);
        m_native_parameter_thread_storage.resize(num_threads);
        m_surface_parameter_thread_storage.resize(num_threads);
        for (auto& [name, gradient] : m_derivative_gradients) {
            gradient.setZero();
        }
        for (auto& [name, gradient] : m_surface_gradients) {
            gradient.setZero();
        }
        for (int thread = 0; thread < num_threads; ++thread) {
            m_native_thread_storage[thread].resize(NSTOKES, this->m_ngeometry);
            m_native_parameter_thread_storage[thread].resize(
                this->m_ngeometry, this->m_wavelength_block_capacity);
            m_surface_parameter_thread_storage[thread].resize(
                this->m_wavelength_block_capacity);
            for (const auto& [name, gradient] : m_derivative_gradients) {
                m_thread_derivative_gradients[thread][name] =
                    Eigen::VectorXd::Zero(gradient.size());
            }
            for (const auto& [name, gradient] : m_surface_gradients) {
                m_thread_surface_gradients[thread][name] =
                    Eigen::VectorXd::Zero(gradient.size());
            }
        }
    }

    template <int NSTOKES>
    Eigen::Vector<double, NSTOKES>
    OutputVJP<NSTOKES>::native_cotangent(int losidx, int wavelidx) const {
        const int linear_index =
            NSTOKES * this->m_nlos * wavelidx + NSTOKES * losidx;
        Eigen::Vector<double, NSTOKES> cotangent;
        for (int stokes = 0; stokes < NSTOKES; ++stokes) {
            cotangent(stokes) = m_cotangent(linear_index + stokes);
        }
        double c = 1.0;
        double s = 0.0;
        if constexpr (NSTOKES >= 3) {
            c = this->m_stokes_C[losidx];
            s = this->m_stokes_S[losidx];
        }
        return transpose_rotate_stokes<NSTOKES>(cotangent, c, s);
    }

    template <int NSTOKES>
    void OutputVJP<NSTOKES>::assign_native_value(
        int losidx, int wavelidx, const Eigen::Vector<double, NSTOKES>& value) {
        const int linear_index =
            NSTOKES * this->m_nlos * wavelidx + NSTOKES * losidx;
        double c = 1.0;
        double s = 0.0;
        if constexpr (NSTOKES >= 3) {
            c = this->m_stokes_C[losidx];
            s = this->m_stokes_S[losidx];
        }
        const auto rotated = rotate_stokes<NSTOKES>(value, c, s);
        for (int stokes = 0; stokes < NSTOKES; ++stokes) {
            m_radiance(linear_index + stokes) = rotated(stokes);
        }
    }

    template <int NSTOKES>
    void OutputVJP<NSTOKES>::accumulate_native_gradient(
        int wavelidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_gradient) {
        Eigen::MatrixXd gradient_block(native_gradient.size(), 1);
        gradient_block.col(0) = native_gradient;
        accumulate_native_gradient(sasktran2::WavelengthBlock<>{wavelidx, 1},
                                   threadidx, gradient_block);
    }

    template <int NSTOKES>
    void OutputVJP<NSTOKES>::accumulate_native_gradient(
        const sasktran2::WavelengthBlock<>& block, int threadidx,
        const Eigen::MatrixXd& native_gradient) {
        if (threadidx < 0 ||
            threadidx >= m_native_parameter_thread_storage.size() ||
            block.count < 1 ||
            block.count > this->m_wavelength_block_capacity ||
            native_gradient.rows() != this->m_atmosphere->num_deriv() ||
            native_gradient.cols() < block.count) {
            throw std::invalid_argument(
                "Invalid native VJP gradient wavelength block");
        }
        const int num_geometry = this->m_ngeometry;
        auto native_parameter =
            m_native_parameter_thread_storage[threadidx].leftCols(block.count);
        for (const auto& [name, gradient] : m_derivative_gradients) {
            const auto& mapping =
                this->m_atmosphere->storage().derivative_mappings_const().at(
                    name);
            const auto& native_mapping = mapping.native_mapping();
            native_parameter.setZero();
            if (native_mapping.d_extinction.has_value()) {
                native_parameter.array() +=
                    native_mapping.d_extinction.value()
                        .middleCols(block.start, block.count)
                        .array() *
                    native_gradient.block(0, 0, num_geometry, block.count)
                        .array();
            }
            if (native_mapping.d_ssa.has_value()) {
                native_parameter.array() +=
                    native_mapping.d_ssa.value()
                        .middleCols(block.start, block.count)
                        .array() *
                    native_gradient
                        .block(this->m_atmosphere->ssa_deriv_start_index(), 0,
                               num_geometry, block.count)
                        .array();
            }
            if (mapping.is_scattering_derivative()) {
                const int start = this->m_atmosphere->scat_deriv_start_index() +
                                  mapping.get_scattering_index() * num_geometry;
                native_parameter.array() +=
                    native_mapping.scat_factor.value()
                        .middleCols(block.start, block.count)
                        .array() *
                    native_gradient.block(start, 0, num_geometry, block.count)
                        .array();
            }
            if (native_mapping.d_emission.has_value() &&
                this->m_atmosphere->include_emission_derivatives()) {
                native_parameter.array() +=
                    native_mapping.d_emission.value()
                        .middleCols(block.start, block.count)
                        .array() *
                    native_gradient
                        .block(this->m_atmosphere->emission_deriv_start_index(),
                               0, num_geometry, block.count)
                        .array();
            }
            auto& target = m_thread_derivative_gradients[threadidx].at(name);
            if (mapping.get_interpolator_const().has_value()) {
                target.noalias() +=
                    mapping.get_interpolator_const().value().transpose() *
                    native_parameter.rowwise().sum();
            } else {
                target += native_parameter.rowwise().sum();
            }
        }

        for (const auto& [name, gradient] : m_surface_gradients) {
            const auto& mapping =
                this->m_atmosphere->surface().derivative_mappings().at(name);
            const auto& native_mapping = mapping.native_surface_mapping();
            auto native_parameter_gradient =
                m_surface_parameter_thread_storage[threadidx].head(block.count);
            native_parameter_gradient.setZero();
            if (native_mapping.d_brdf.has_value()) {
                const int num_surface_derivatives =
                    this->m_atmosphere->surface().num_deriv();
                for (int lane = 0; lane < block.count; ++lane) {
                    native_parameter_gradient(lane) +=
                        native_mapping.d_brdf.value()
                            .row(block.wavelength(lane))
                            .dot(native_gradient.col(lane).segment(
                                this->m_atmosphere->surface_deriv_start_index(),
                                num_surface_derivatives));
                }
            }
            if (native_mapping.d_emission.has_value() &&
                this->m_atmosphere->include_emission_derivatives()) {
                for (int lane = 0; lane < block.count; ++lane) {
                    native_parameter_gradient(lane) +=
                        native_mapping.d_emission.value()(
                            block.wavelength(lane), 0) *
                        native_gradient(
                            this->m_atmosphere
                                ->surface_emission_deriv_start_index(),
                            lane);
                }
            }
            auto& target = m_thread_surface_gradients[threadidx].at(name);
            if (mapping.get_interpolator_const().has_value()) {
                target.noalias() += mapping.get_interpolator_const()
                                        .value()
                                        .middleRows(block.start, block.count)
                                        .transpose() *
                                    native_parameter_gradient;
            } else {
                target.segment(block.start, block.count) +=
                    native_parameter_gradient.transpose();
            }
        }
    }

    template <int NSTOKES>
    template <typename Radiance>
    void OutputVJP<NSTOKES>::assign_lane(const Radiance& radiance, int losidx,
                                         int wavelidx, int threadidx) {
        const int linear_index =
            NSTOKES * this->m_nlos * wavelidx + NSTOKES * losidx;
        double c = 1.0;
        double s = 0.0;
        if constexpr (NSTOKES >= 3) {
            c = this->m_stokes_C[losidx];
            s = this->m_stokes_S[losidx];
        }

        Eigen::Matrix<double, NSTOKES, 1> value;
        Eigen::Matrix<double, NSTOKES, 1> external_cotangent;
        for (int stokes = 0; stokes < NSTOKES; ++stokes) {
            value(stokes) = radiance.value(stokes);
            external_cotangent(stokes) = m_cotangent(linear_index + stokes);
        }
        const auto internal_cotangent =
            transpose_rotate_stokes<NSTOKES>(external_cotangent, c, s);
        value = rotate_stokes<NSTOKES>(value, c, s);
        for (int stokes = 0; stokes < NSTOKES; ++stokes) {
            m_radiance(linear_index + stokes) = value(stokes);
        }

        auto& native = m_native_thread_storage[threadidx];
        for (const auto& [name, gradient] : m_derivative_gradients) {
            const auto& mapping =
                this->m_atmosphere->storage().derivative_mappings_const().at(
                    name);
            mapped_volume_derivatives<NSTOKES>(
                radiance, mapping, wavelidx, this->m_ngeometry,
                this->m_atmosphere->num_scattering_deriv_groups(),
                this->m_atmosphere->include_emission_derivatives(), native);
            const Eigen::VectorXd native_gradient =
                native.transpose() * internal_cotangent;
            auto& thread_gradient =
                m_thread_derivative_gradients[threadidx].at(name);
            if (mapping.get_interpolator_const().has_value()) {
                const auto& interpolator =
                    mapping.get_interpolator_const().value();
                if (thread_gradient.size() != interpolator.cols()) {
                    throw std::invalid_argument(
                        "Volume linearization gradient has incorrect size");
                }
                thread_gradient.noalias() +=
                    interpolator.transpose() * native_gradient;
            } else {
                if (thread_gradient.size() != native_gradient.size()) {
                    throw std::invalid_argument(
                        "Volume linearization gradient has incorrect size");
                }
                thread_gradient += native_gradient;
            }
        }

        for (const auto& [name, gradient] : m_surface_gradients) {
            const auto& mapping =
                this->m_atmosphere->surface().derivative_mappings().at(name);
            const auto native_derivative = mapped_surface_derivative<NSTOKES>(
                radiance, mapping, wavelidx, this->m_ngeometry,
                this->m_atmosphere->num_source_deriv_groups(),
                this->m_atmosphere->surface().num_deriv(),
                this->m_atmosphere->surface_emission_deriv_start_index(),
                this->m_atmosphere->include_emission_derivatives());
            const double native_gradient =
                native_derivative.dot(internal_cotangent);
            auto& thread_gradient =
                m_thread_surface_gradients[threadidx].at(name);
            if (mapping.get_interpolator_const().has_value()) {
                const auto& interpolator =
                    mapping.get_interpolator_const().value();
                if (thread_gradient.size() != interpolator.cols()) {
                    throw std::invalid_argument(
                        "Surface linearization gradient has incorrect size");
                }
                thread_gradient.noalias() +=
                    interpolator.row(wavelidx).transpose() * native_gradient;
            } else {
                if (thread_gradient.size() != this->m_nwavel) {
                    throw std::invalid_argument(
                        "Native spectral surface gradient has incorrect size");
                }
                thread_gradient(wavelidx) += native_gradient;
            }
        }
    }

    template <int NSTOKES>
    void OutputVJP<NSTOKES>::assign(
        const sasktran2::WavelengthBlock<>& block,
        const sasktran2::WavelengthBlockDual<NSTOKES>& radiance, int losidx,
        int threadidx) {
        for (int lane = 0; lane < block.count; ++lane) {
            const sasktran2::WavelengthBlockConstLaneDualView<NSTOKES>
                radiance_lane(radiance, lane);
            assign_lane(radiance_lane, losidx, block.start + lane, threadidx);
        }
    }

    template <int NSTOKES> void OutputVJP<NSTOKES>::finalize() {
        for (auto& [name, gradient] : m_derivative_gradients) {
            gradient.setZero();
            for (const auto& thread_gradients : m_thread_derivative_gradients) {
                gradient += thread_gradients.at(name);
            }
        }
        for (auto& [name, gradient] : m_surface_gradients) {
            gradient.setZero();
            for (const auto& thread_gradients : m_thread_surface_gradients) {
                gradient += thread_gradients.at(name);
            }
        }
    }

    template class OutputJVP<1>;
    template class OutputJVP<3>;
    template class OutputVJP<1>;
    template class OutputVJP<3>;
} // namespace sasktran2
