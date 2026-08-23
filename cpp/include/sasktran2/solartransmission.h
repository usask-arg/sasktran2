#pragma once

#include "sasktran2/atmosphere/grid_storage.h"
#include "sasktran2/geometry.h"
#include <sasktran2/source_interface.h>
#include <sasktran2/raytracing.h>
#include <sasktran2/atmosphere/atmosphere.h>
#include <sasktran2/config.h>
#include <sasktran2/dual.h>
#include <array>
#include <cstdint>
#include <limits>

namespace sasktran2::solartransmission {
    /** Generates one row of what is known as the geometry matrix.  The optical
     * depth at every necessary integration point can be written as \f$M k\f$
     * where \f$k\f$ is the extinction coefficient vector.  In 1D geometry, the
     *  matrix \f$M\f$ is densish, but it may be sparse in higher dimensions.
     * This function constructs one row of the matrix \f$M\f$, which corresponds
     * to one traced ray.
     *
     * @param row row index
     * @param traced_ray ray traced to the sun
     * @param result matrix to store the result in
     */
    inline void assign_dense_matrix_column(
        int row, const sasktran2::raytracing::TracedRay& traced_ray,
        Eigen::MatrixXd& result) {
        for (int i = 0; i < traced_ray.layers.size(); ++i) {
            const auto weights = traced_ray.optical_depth_weights(i);
            for (std::size_t j = 0; j < weights.size(); ++j) {
                const auto weight = weights[j];
                result(row, weight.first) += weight.second;
            }
        }
    }

    /** Row-major solar geometry storage with a compact exact-2D mode.
     *
     * Exact 2D calculations commonly have hundreds of thousands of rows but
     * only a few thousand atmosphere columns. Eigen uses the same index type
     * for both dimensions, so its inner indices remain 32 bit. The compact
     * mode keeps row offsets as size_t while storing column indices as uint16.
     * Grids wider than the compact index range transparently use Eigen's
     * standard sparse representation.
     */
    class SolarGeometryMatrix {
      public:
        using StandardMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;

        class InnerIterator {
          private:
            const double* m_values = nullptr;
            const std::uint16_t* m_compact_indices = nullptr;
            const int* m_standard_indices = nullptr;
            std::size_t m_position = 0;
            std::size_t m_end = 0;

          public:
            InnerIterator(const SolarGeometryMatrix& matrix, Eigen::Index row) {
                if (matrix.m_compact) {
                    m_position = matrix.m_compact_outer[row];
                    m_end = matrix.m_compact_outer[row + 1];
                    m_values = matrix.m_compact_values.data();
                    m_compact_indices = matrix.m_compact_inner.data();
                } else {
                    const auto* outer = matrix.m_standard.outerIndexPtr();
                    m_position = static_cast<std::size_t>(outer[row]);
                    const auto* inner_nonzeros =
                        matrix.m_standard.innerNonZeroPtr();
                    m_end = inner_nonzeros == nullptr
                                ? static_cast<std::size_t>(outer[row + 1])
                                : m_position + static_cast<std::size_t>(
                                                   inner_nonzeros[row]);
                    m_values = matrix.m_standard.valuePtr();
                    m_standard_indices = matrix.m_standard.innerIndexPtr();
                }
            }

            explicit operator bool() const { return m_position < m_end; }
            InnerIterator& operator++() {
                ++m_position;
                return *this;
            }
            Eigen::Index index() const {
                return m_compact_indices == nullptr
                           ? m_standard_indices[m_position]
                           : m_compact_indices[m_position];
            }
            double value() const { return m_values[m_position]; }
        };

      private:
        StandardMatrix m_standard;
        std::vector<double> m_compact_values;
        std::vector<std::uint16_t> m_compact_inner;
        std::vector<std::size_t> m_compact_outer;
        Eigen::Index m_compact_rows = 0;
        Eigen::Index m_compact_cols = 0;
        bool m_compact = false;

        void reserve_capacity(Eigen::Index capacity) {
            if (m_compact) {
                m_compact_values.reserve(static_cast<std::size_t>(capacity));
                m_compact_inner.reserve(static_cast<std::size_t>(capacity));
            } else {
                auto& storage = m_standard.data();
                storage.reserve(capacity - storage.size());
            }
        }

      public:
        StandardMatrix& use_standard() {
            m_compact = false;
            std::vector<double>().swap(m_compact_values);
            std::vector<std::uint16_t>().swap(m_compact_inner);
            std::vector<std::size_t>().swap(m_compact_outer);
            m_compact_rows = 0;
            m_compact_cols = 0;
            return m_standard;
        }

        const StandardMatrix& standard() const { return m_standard; }

        void initialize_exact(Eigen::Index rows, Eigen::Index cols,
                              Eigen::Index initial_capacity) {
            m_compact = cols <= static_cast<Eigen::Index>(
                                    std::numeric_limits<std::uint16_t>::max()) +
                                    1;
            if (!m_compact) {
                auto& standard = use_standard();
                standard.resize(rows, cols);
                standard.reserve(initial_capacity);
                return;
            }

            m_standard = StandardMatrix{};
            m_compact_rows = rows;
            m_compact_cols = cols;
            m_compact_values.clear();
            m_compact_inner.clear();
            m_compact_outer.assign(static_cast<std::size_t>(rows) + 1, 0);
            reserve_capacity(initial_capacity);
        }

        void ensure_capacity(Eigen::Index additional) {
            const Eigen::Index required = non_zeros() + additional;
            const Eigen::Index allocated = allocated_size();
            if (required <= allocated) {
                return;
            }
            const Eigen::Index grown = allocated + allocated / 2;
            reserve_capacity(std::max(required, grown));
        }

        void start_row(Eigen::Index row) {
            if (m_compact) {
                m_compact_outer[row] = m_compact_values.size();
            } else {
                m_standard.startVec(row);
            }
        }

        void insert_back(Eigen::Index row, Eigen::Index column, double value) {
            if (m_compact) {
                m_compact_inner.push_back(static_cast<std::uint16_t>(column));
                m_compact_values.push_back(value);
            } else {
                m_standard.insertBackByOuterInner(row, column) = value;
            }
        }

        void finalize() {
            if (m_compact) {
                m_compact_outer[m_compact_rows] = m_compact_values.size();
            } else {
                m_standard.finalize();
                m_standard.data().squeeze();
            }
        }

        Eigen::Index rows() const {
            return m_compact ? m_compact_rows : m_standard.rows();
        }
        Eigen::Index cols() const {
            return m_compact ? m_compact_cols : m_standard.cols();
        }
        Eigen::Index non_zeros() const {
            return m_compact
                       ? static_cast<Eigen::Index>(m_compact_values.size())
                       : m_standard.nonZeros();
        }
        Eigen::Index allocated_size() const {
            return m_compact
                       ? static_cast<Eigen::Index>(m_compact_values.capacity())
                       : m_standard.data().allocatedSize();
        }
        Eigen::Index row_nonzeros(Eigen::Index row) const {
            if (m_compact) {
                return static_cast<Eigen::Index>(m_compact_outer[row + 1] -
                                                 m_compact_outer[row]);
            }
            return m_standard.innerVector(row).nonZeros();
        }
        bool is_compact() const { return m_compact; }

        template <typename Vector>
        double row_dot(Eigen::Index row, const Vector& vector) const {
            double result = 0.0;
            if (m_compact) {
                const auto begin = m_compact_outer[row];
                const auto end = m_compact_outer[row + 1];
                for (std::size_t position = begin; position < end; ++position) {
                    result += m_compact_values[position] *
                              vector(m_compact_inner[position]);
                }
                return result;
            }

            for (StandardMatrix::InnerIterator entry(m_standard, row); entry;
                 ++entry) {
                result += entry.value() * vector(entry.index());
            }
            return result;
        }

        template <typename Vector>
        void accumulate_row(Eigen::Index row, double scale,
                            Vector& vector) const {
            if (m_compact) {
                const auto begin = m_compact_outer[row];
                const auto end = m_compact_outer[row + 1];
                for (std::size_t position = begin; position < end; ++position) {
                    vector(m_compact_inner[position]) +=
                        scale * m_compact_values[position];
                }
                return;
            }

            for (StandardMatrix::InnerIterator entry(m_standard, row); entry;
                 ++entry) {
                vector(entry.index()) += scale * entry.value();
            }
        }

        template <typename Rhs, typename Destination>
        void multiply(const Rhs& rhs, Destination& destination) const {
            if (!m_compact) {
                destination.noalias() = m_standard * rhs;
                return;
            }

            if (rhs.cols() == 1) {
                for (Eigen::Index row = 0; row < m_compact_rows; ++row) {
                    double result = 0.0;
                    const auto begin = m_compact_outer[row];
                    const auto end = m_compact_outer[row + 1];
                    for (std::size_t position = begin; position < end;
                         ++position) {
                        result += m_compact_values[position] *
                                  rhs(m_compact_inner[position], 0);
                    }
                    destination(row, 0) = result;
                }
                return;
            }

            destination.setZero();
            for (Eigen::Index row = 0; row < m_compact_rows; ++row) {
                const auto begin = m_compact_outer[row];
                const auto end = m_compact_outer[row + 1];
                for (std::size_t position = begin; position < end; ++position) {
                    const auto column = m_compact_inner[position];
                    const double value = m_compact_values[position];
                    for (Eigen::Index lane = 0; lane < rhs.cols(); ++lane) {
                        destination(row, lane) += value * rhs(column, lane);
                    }
                }
            }
        }
    };

    /** Compact row-major interpolation from solar-table nodes to ray
     * endpoints.
     *
     * The rows are assembled in order and normally contain at most eight
     * entries for a three-dimensional table. Keeping the construction in CSR
     * form avoids the large temporary triplet list required by Eigen's sparse
     * matrix builder for the hundreds of thousands of endpoints used by the
     * successive-orders source.
     */
    class SolarTableInterpolation {
      private:
        std::vector<std::uint32_t> m_outer;
        std::vector<std::uint32_t> m_inner;
        std::vector<double> m_values;
        Eigen::Index m_rows = 0;
        Eigen::Index m_cols = 0;
        Eigen::Index m_next_row = 0;

      public:
        void clear();
        void initialize(Eigen::Index rows, Eigen::Index cols,
                        Eigen::Index capacity);
        void append_row(
            const std::vector<std::pair<int, double>>& interpolation_weights);
        void finalize();

        Eigen::Index rows() const { return m_rows; }
        Eigen::Index cols() const { return m_cols; }
        Eigen::Index non_zeros() const {
            return static_cast<Eigen::Index>(m_values.size());
        }

        void apply(Eigen::Ref<const Eigen::VectorXd> table_values,
                   Eigen::Ref<Eigen::VectorXd> endpoint_values) const;
        void apply_transpose(Eigen::Ref<const Eigen::VectorXd> endpoint_values,
                             Eigen::Ref<Eigen::VectorXd> table_values) const;
        std::size_t storage_bytes() const;
    };

    /** Geometry-only solar optical-depth table shared by source terms.
     *
     * Implementations may use a conventional geometry matrix or a
     * characteristic sweep. The latter stores only incremental path stencils,
     * which is substantially smaller than materializing cumulative
     * atmosphere-to-sun rows at every source-ray endpoint.
     */
    class SolarTransmissionTableEvaluator {
      public:
        virtual ~SolarTransmissionTableEvaluator() = default;

        virtual void initialize_config(const sasktran2::Config& config) = 0;
        virtual void
        initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>&
                                integration_rays) = 0;
        virtual void generate_interpolation(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            SolarTableInterpolation& interpolator,
            std::vector<bool>& ground_hit_flag,
            std::vector<Eigen::Vector3d>* solar_propagation_directions =
                nullptr) const = 0;

        virtual Eigen::Index table_size() const = 0;
        virtual Eigen::Index atmosphere_size() const = 0;
        virtual void apply(Eigen::Ref<const Eigen::VectorXd> extinction,
                           Eigen::Ref<Eigen::VectorXd> table_od) const = 0;
        virtual void
        accumulate_transpose(Eigen::Ref<const Eigen::VectorXd> table_cotangent,
                             Eigen::Ref<Eigen::VectorXd> extinction_cotangent,
                             double scale) const = 0;
        virtual std::size_t storage_bytes() const = 0;
    };

    class SolarTransmissionBase {
      protected:
        const Geometry& m_geometry;
        const Geometry1D* m_geometry_1d = nullptr;
        const sasktran2::raytracing::RayTracerBase* m_raytracer = nullptr;
#ifdef SKTRAN_RUST_SUPPORT
        const Geometry2D* m_geometry_2d = nullptr;
        const sasktran2::raytracing::RustRayTracer2D* m_raytracer_2d = nullptr;
#endif

      public:
        SolarTransmissionBase(
            const Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer)
            : m_geometry(geometry), m_geometry_1d(&geometry),
              m_raytracer(&raytracer) {}

#ifdef SKTRAN_RUST_SUPPORT
        SolarTransmissionBase(
            const Geometry2D& geometry,
            const sasktran2::raytracing::RustRayTracer2D& raytracer)
            : m_geometry(geometry), m_geometry_2d(&geometry),
              m_raytracer_2d(&raytracer) {}
#endif
    };

    class SolarTransmissionExact : public SolarTransmissionBase {
      private:
      public:
        SolarTransmissionExact(
            const Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer)
            : SolarTransmissionBase(geometry, raytracer) {}

        virtual void initialize_config(const sasktran2::Config& config){};

        virtual void
        initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>&
                                integration_rays){};

        void generate_geometry_matrix(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            Eigen::MatrixXd& od_matrix,
            std::vector<bool>& ground_hit_flag) const;

#ifdef SKTRAN_RUST_SUPPORT
        SolarTransmissionExact(
            const Geometry2D& geometry,
            const sasktran2::raytracing::RustRayTracer2D& raytracer)
            : SolarTransmissionBase(geometry, raytracer) {}

        void generate_geometry_matrix(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            SolarGeometryMatrix& od_matrix,
            std::vector<bool>& ground_hit_flag) const;

        void generate_refracted_geometry_matrix(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            const std::vector<Eigen::Vector3d>& solar_propagation_directions,
            SolarGeometryMatrix& od_matrix,
            std::vector<bool>& ground_hit_flag) const;
#endif
    };

    class SolarTransmissionTable : public SolarTransmissionExact,
                                   public SolarTransmissionTableEvaluator {
      private:
        std::unique_ptr<sasktran2::grids::SourceLocationInterpolator>
            m_location_interpolator;
        const sasktran2::Config* m_config;
        Eigen::MatrixXd m_geometry_matrix;

        std::vector<bool> m_ground_hit_flag;

      public:
        SolarTransmissionTable(
            const Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer)
            : SolarTransmissionExact(geometry, raytracer) {}

        void initialize_config(const sasktran2::Config& config) override {
            m_config = &config;
        };

        void
        initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>&
                                integration_rays) override;

        void generate_interpolation_matrix(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            Eigen::SparseMatrix<double, Eigen::RowMajor>& interpolator,
            std::vector<bool>& ground_hit_flag) const;

        void generate_interpolation(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            SolarTableInterpolation& interpolator,
            std::vector<bool>& ground_hit_flag,
            std::vector<Eigen::Vector3d>* solar_propagation_directions =
                nullptr) const override;
        Eigen::Index table_size() const override {
            return m_geometry_matrix.rows();
        }
        Eigen::Index atmosphere_size() const override {
            return m_geometry_matrix.cols();
        }
        void apply(Eigen::Ref<const Eigen::VectorXd> extinction,
                   Eigen::Ref<Eigen::VectorXd> table_od) const override;
        void
        accumulate_transpose(Eigen::Ref<const Eigen::VectorXd> table_cotangent,
                             Eigen::Ref<Eigen::VectorXd> extinction_cotangent,
                             double scale) const override;
        std::size_t storage_bytes() const override;

        const Eigen::MatrixXd& geometry_matrix() const {
            return m_geometry_matrix;
        }
    };

#ifdef SKTRAN_RUST_SUPPORT
    /** Three-dimensional solar optical-depth table for a structured 2D
     * atmosphere.
     *
     * Table locations are parameterized by altitude, solar zenith angle, and
     * azimuth around the solar axis. Parallel rays are launched at the
     * top-of-atmosphere chord and traced inward. Each ray stores only the
     * incremental 2D cell-basis stencils; optical depth is accumulated by a
     * forward sweep and its VJP by the corresponding reverse sweep.
     */
    class SolarTransmissionTable2D final
        : public SolarTransmissionTableEvaluator {
      private:
        class Impl;
        std::unique_ptr<Impl> m_impl;

      public:
        SolarTransmissionTable2D(
            const Geometry2D& geometry,
            const sasktran2::raytracing::RustRayTracer2D& raytracer);
        ~SolarTransmissionTable2D() override;

        void initialize_config(const sasktran2::Config& config) override;
        void
        initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>&
                                integration_rays) override;
        void generate_interpolation(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            SolarTableInterpolation& interpolator,
            std::vector<bool>& ground_hit_flag,
            std::vector<Eigen::Vector3d>* solar_propagation_directions =
                nullptr) const override;
        void generate_solar_geometry(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            std::vector<bool>& ground_hit_flag,
            std::vector<Eigen::Vector3d>& solar_propagation_directions) const;
        Eigen::Index table_size() const override;
        Eigen::Index atmosphere_size() const override;
        void apply(Eigen::Ref<const Eigen::VectorXd> extinction,
                   Eigen::Ref<Eigen::VectorXd> table_od) const override;
        void
        accumulate_transpose(Eigen::Ref<const Eigen::VectorXd> table_cotangent,
                             Eigen::Ref<Eigen::VectorXd> extinction_cotangent,
                             double scale) const override;
        std::size_t storage_bytes() const override;
    };
#endif

    /**
     * The PhaseHandler is responsible for constructing the phase function for
     * the single scatter source term. The single scatter source term is needed
     * on a set of (wavelength, cos_angle, geometry) points.
     *
     * Construction can be done in one of two ways, depending on the
     * configuration and user input. The default is to construct the phase
     * function through the user input Legendre coefficients.
     *
     * The other option is for the user to directly input the phase function.
     * This is useful for cases where the phase function requires many terms in
     * the Legendre series.
     *
     * Storage for the phase function is handled in a semi-complicated hard to
     * understand manner. For each thread, we store the phase function (stokes,
     * internal_index) where internal_index represents a single scattering angle
     * at a single geomtry level in the atmosphere.
     *
     * The internal_index is mapped back to the geometry index through the
     * m_internal_to_geometry, and the scattering angle is determined by the
     * m_internal_to_cos_scatter.  This combination lets us calculate the phase
     * function.  But to actually use it, we need to map the entrance and exit
     * points of the ray to the internal indices.  This is done through the
     * packed entrance/exit ranges and m_geometry_to_internal
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES> class PhaseHandler {
      private:
        using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic,
                                             Eigen::Dynamic, Eigen::RowMajor>;

        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere;
        const sasktran2::Config* m_config;
        const sasktran2::Geometry& m_geometry;

        std::vector<std::array<double, NSTOKES>>
            m_scatter_angles; /** Full list of  scattering angles that we need.
                                 for NSTOKES =3 this is (cos_scatter, C1, C2) */

        Eigen::MatrixXd m_wigner_d00; /** Wigner D matrix for the phase function
                                         (legendre_order, scatter_angle) */
        Eigen::MatrixXd m_wigner_d02; /** Wigner D matrix for the phase function
                                         (legendre_order, scatter_angle) */

        // Internal phase functions and derivatives on the actual grid
        Eigen::Tensor<double, 3>
            m_phase; /** (stokes eq, internal_index, thread) **/
        Eigen::Tensor<double, 4>
            m_d_phase; /** (stokes eq, internal_index, deriv, thread) **/

        // Batch storage is wavelength-contiguous. The row mappings are
        // phase_component * num_internal + internal_index and
        // (deriv * num_phase_components + phase_component) * num_internal +
        // internal_index, respectively.
        std::vector<RowMajorMatrix> m_phase_batch;
        std::vector<RowMajorMatrix> m_d_phase_batch;
        int m_wavelength_batch_capacity = 0;
        std::vector<int> m_active_wavelength_block_start;
        std::vector<int> m_active_wavelength_block_count;

        std::vector<std::uint32_t> m_geometry_layer_offsets;
        std::vector<std::uint32_t> m_geometry_entrance_offsets;
        std::vector<std::uint32_t> m_geometry_exit_offsets;
        std::vector<int> m_geometry_to_internal;
        std::vector<int> m_internal_to_geometry; /** Maps the internal index
                                                      to the geometry index */
        std::vector<int>
            m_internal_to_cos_scatter; /** Determines what scattering angle to
                                          use for each internal index */

      public:
        PhaseHandler(const Geometry& geometry) : m_geometry(geometry) {}

        /**
         * Initializes the phase handler with the configuration object
         */
        void initialize_config(const sasktran2::Config& config) {
            m_config = &config;
        }

        /**
         *  Initializes the phase handler with the atmosphere object
         *
         *  @param atmosphere The atmosphere object
         */
        void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere);

        /**
         * Initializes the phase handler with the geometry object
         *
         * @param los_rays The traced line of sight rays
         * @param index_map The index map
         */
        void initialize_geometry(
            const std::vector<sasktran2::raytracing::TracedRay>& los_rays,
            const std::vector<std::vector<int>>& index_map,
            const std::vector<Eigen::Vector3d>* solar_propagation_directions =
                nullptr);

        /**
         *   Calculates the phase function from the legendre coefficients at the
         * necessary scatter angles
         *
         *   @param threadidx The thread index
         *   @param wavelidx The wavelength index
         */
        void calculate(int wavelidx, int threadidx);

        void initialize_wavelength_blocks(int batch_size);

        template <int N>
        void calculate_block(const sasktran2::WavelengthBlock<N>& batch,
                             int threadidx);

        /**
         * Calculates the phase function at a given point and puts it into
         * source
         *
         * @param wavelidx The wavelength index
         * @param losidx The line of sight index
         * @param layeridx The layer index
         * @param index_weights The interpolation weights
         * @param is_entrance  True if we are at the entrance to a layer, false
         * if we are at the exit to a layer
         * @param source The source term
         */
        void scatter(int wavelidx, int losidx, int layeridx,
                     const raytracing::GridWeightStencilView& index_weights,
                     bool is_entrance,
                     sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                     NSTOKES>& source) const;

        /** Returns the interpolated phase vector and accumulates its
         * scattering-property derivatives directly into a target source. */
        Eigen::Vector<double, NSTOKES> scatter_and_accumulate_derivative(
            int threadidx, int losidx, int layeridx,
            const raytracing::GridWeightStencilView& index_weights,
            bool is_entrance, double source_amplitude, double derivative_scale,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                target) const;

        Eigen::Vector<double, NSTOKES> scatter_and_accumulate_derivative(
            int threadidx, int losidx, int layeridx,
            const raytracing::GridWeightStencilView& index_weights,
            bool is_entrance, double source_amplitude, double derivative_scale,
            sasktran2::WavelengthBlockLaneDualView<NSTOKES, 1>& target) const;

        template <int N>
        void scatter_and_accumulate_derivative_block(
            int threadidx, int losidx, int layeridx,
            const raytracing::GridWeightStencilView& index_weights,
            bool is_entrance, const sasktran2::WavelengthBlock<N>& batch,
            const Eigen::Ref<
                const Eigen::Matrix<double, 1, N, Eigen::RowMajor>>&
                source_amplitude,
            const Eigen::Ref<
                const Eigen::Matrix<double, 1, N, Eigen::RowMajor>>&
                derivative_scale,
            sasktran2::WavelengthBlockDual<NSTOKES>& target,
            Eigen::Matrix<double, NSTOKES, Eigen::Dynamic, Eigen::RowMajor>&
                phase_result) const;

        void scatter_jvp(int threadidx, int losidx, int layeridx, int wavelidx,
                         const raytracing::GridWeightStencilView& index_weights,
                         bool is_entrance,
                         Eigen::Ref<const Eigen::VectorXd> native_tangent,
                         Eigen::Vector<double, NSTOKES>& phase,
                         Eigen::Vector<double, NSTOKES>& phase_jvp) const;

        Eigen::Vector<double, NSTOKES>
        scatter_value(int threadidx, int losidx, int layeridx, int wavelidx,
                      const raytracing::GridWeightStencilView& index_weights,
                      bool is_entrance) const;

        void scatter_vjp(int threadidx, int losidx, int layeridx, int wavelidx,
                         const raytracing::GridWeightStencilView& index_weights,
                         bool is_entrance,
                         const Eigen::Vector<double, NSTOKES>& phase_cotangent,
                         Eigen::Ref<Eigen::VectorXd> native_gradient) const;

      private:
        template <typename Target>
        Eigen::Vector<double, NSTOKES> scatter_and_accumulate_derivative_impl(
            int threadidx, int losidx, int layeridx,
            const raytracing::GridWeightStencilView& index_weights,
            bool is_entrance, double source_amplitude, double derivative_scale,
            Target& target) const;

        void initialize_geometry_impl(
            const std::vector<sasktran2::raytracing::TracedRay>& los_rays,
            const std::vector<std::vector<int>>& index_map,
            const std::vector<Eigen::Vector3d>* solar_propagation_directions);

        const int* geometry_internal_indices(int losidx, int layeridx,
                                             bool is_entrance) const {
            const auto flat_layer = m_geometry_layer_offsets[losidx] +
                                    static_cast<std::uint32_t>(layeridx);
            const auto offset = is_entrance
                                    ? m_geometry_entrance_offsets[flat_layer]
                                    : m_geometry_exit_offsets[flat_layer];
            return m_geometry_to_internal.data() + offset;
        }

        void
        scatter_impl(int wavelidx, int losidx, int layeridx,
                     const raytracing::GridWeightStencilView& index_weights,
                     bool is_entrance,
                     sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                     NSTOKES>& source) const;
    };

    template <int NSTOKES>
    inline void
    scattering_source(const PhaseHandler<NSTOKES>& phase_handler, int threadidx,
                      int losidx, int layeridx, int wavelidx,
                      const raytracing::GridWeightStencilView& index_weights,
                      bool is_entrance, double solar_trans,
                      const atmosphere::Atmosphere<NSTOKES>& atmosphere,
                      SolarGeometryMatrix::InnerIterator solar_trans_iter,
                      bool calculate_derivatives,
                      sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                      NSTOKES>& source) {
        const auto& storage = atmosphere.storage();
        source.value.setZero();

        if (calculate_derivatives) {
            source.deriv.setZero();
        }

        double ssa = 0;
        double k = 0;
        if (index_weights.size() == 2) {
            // Structured 1D layers always have two backing nodes. Keep the
            // common ray representation while retaining the direct endpoint
            // interpolation used by the previous 1D-specific stencil.
            const auto lower = index_weights[0];
            const auto upper = index_weights[1];
            if (lower.second == 0.0) {
                ssa = storage.ssa(upper.first, wavelidx) * upper.second;
                k = storage.total_extinction(upper.first, wavelidx) *
                    upper.second;
            } else if (upper.second == 0.0) {
                ssa = storage.ssa(lower.first, wavelidx) * lower.second;
                k = storage.total_extinction(lower.first, wavelidx) *
                    lower.second;
            } else {
                ssa = storage.ssa(lower.first, wavelidx) * lower.second +
                      storage.ssa(upper.first, wavelidx) * upper.second;
                k = storage.total_extinction(lower.first, wavelidx) *
                        lower.second +
                    storage.total_extinction(upper.first, wavelidx) *
                        upper.second;
            }
        } else {
            for (std::size_t index = 0; index < index_weights.size(); ++index) {
                const auto ele = index_weights[index];
                if (ele.second == 0.0) {
                    continue;
                }
                ssa += storage.ssa(ele.first, wavelidx) * ele.second;
                k += storage.total_extinction(ele.first, wavelidx) * ele.second;
            }
        }

        const double source_amplitude = k * ssa * solar_trans / (EIGEN_PI * 4);
        const bool use_zero_safe_derivative_path =
            calculate_derivatives && (k == 0.0 || ssa == 0.0);
        if (!use_zero_safe_derivative_path) {
            // The common path lets the phase handler work directly with the
            // scaled source, avoiding a second dense derivative pass.
            source.value(0) = source_amplitude;
            phase_handler.scatter(threadidx, losidx, layeridx, index_weights,
                                  is_entrance, source);

            if (!calculate_derivatives) {
                return;
            }
            for (auto it = solar_trans_iter; it; ++it) {
                source.deriv(Eigen::placeholders::all, it.index()) -=
                    it.value() * source.value;
            }
            for (std::size_t index = 0; index < index_weights.size(); ++index) {
                const auto ele = index_weights[index];
                if (ele.second == 0.0) {
                    continue;
                }
                source.deriv(Eigen::placeholders::all,
                             atmosphere.ssa_deriv_start_index() + ele.first) +=
                    ele.second * source.value / ssa;
                source.deriv(Eigen::placeholders::all, ele.first) +=
                    ele.second * source.value / k;
            }
            return;
        }

        // For derivatives at k == 0 or SSA == 0, evaluate the unit-amplitude
        // phase first so the nonzero boundary derivative remains well-defined.
        source.value(0) = 1.0;

        phase_handler.scatter(threadidx, losidx, layeridx, index_weights,
                              is_entrance, source);

        const Eigen::Vector<double, NSTOKES> phase = source.value;
        source.value *= source_amplitude;

        if (!calculate_derivatives) {
            return;
        }
        source.deriv *= source_amplitude;
        // Solar transmission derivative factors
        for (auto it = solar_trans_iter; it; ++it) {
            source.deriv(Eigen::placeholders::all, it.index()) -=
                it.value() * source.value;
        }

        // And SSA/k derivative factors
        for (std::size_t index = 0; index < index_weights.size(); ++index) {
            const auto ele = index_weights[index];
            if (ele.second == 0.0) {
                continue;
            }
            source.deriv(Eigen::placeholders::all,
                         atmosphere.ssa_deriv_start_index() + ele.first) +=
                ele.second * k * solar_trans / (EIGEN_PI * 4) * phase;

            source.deriv(Eigen::placeholders::all, ele.first) +=
                ele.second * ssa * solar_trans / (EIGEN_PI * 4) * phase;
        }
    }

    /** Evaluates one exact single-scatter endpoint and accumulates its
     * derivatives directly into the integrated ray source. This avoids
     * materializing and then copying a dense endpoint derivative buffer. */
    template <int NSTOKES, typename Target>
    inline Eigen::Vector<double, NSTOKES> accumulate_exact_scattering_source(
        const PhaseHandler<NSTOKES>& phase_handler, int threadidx, int losidx,
        int layeridx, int wavelidx,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance, double solar_trans,
        const atmosphere::Atmosphere<NSTOKES>& atmosphere,
        SolarGeometryMatrix::InnerIterator solar_trans_iter,
        double derivative_scale, Target& target) {
        const auto& storage = atmosphere.storage();
        double ssa = 0.0;
        double extinction = 0.0;
        for (std::size_t index = 0; index < index_weights.size(); ++index) {
            const auto weight = index_weights[index];
            if (weight.second == 0.0) {
                continue;
            }
            ssa += storage.ssa(weight.first, wavelidx) * weight.second;
            extinction += storage.total_extinction(weight.first, wavelidx) *
                          weight.second;
        }

        const double unscaled_amplitude = solar_trans / (EIGEN_PI * 4);
        const double source_amplitude = extinction * ssa * unscaled_amplitude;
        const Eigen::Vector<double, NSTOKES> phase =
            phase_handler.scatter_and_accumulate_derivative(
                threadidx, losidx, layeridx, index_weights, is_entrance,
                source_amplitude, derivative_scale, target);
        const Eigen::Vector<double, NSTOKES> endpoint_source =
            source_amplitude * phase;

        for (auto it = solar_trans_iter; it; ++it) {
            target.deriv.col(it.index()) -=
                derivative_scale * it.value() * endpoint_source;
        }

        for (std::size_t index = 0; index < index_weights.size(); ++index) {
            const auto weight = index_weights[index];
            if (weight.second == 0.0) {
                continue;
            }
            target.deriv.col(atmosphere.ssa_deriv_start_index() +
                             weight.first) += derivative_scale * weight.second *
                                              extinction * unscaled_amplitude *
                                              phase;
            target.deriv.col(weight.first) += derivative_scale * weight.second *
                                              ssa * unscaled_amplitude * phase;
        }

        return endpoint_source;
    }

    template <int NSTOKES> struct ExactScatteringBlockScratch {
        using BatchMatrix =
            Eigen::Matrix<double, NSTOKES, Eigen::Dynamic, Eigen::RowMajor>;

        Eigen::RowVectorXd ssa;
        Eigen::RowVectorXd extinction;
        Eigen::RowVectorXd unscaled_amplitude;
        Eigen::RowVectorXd source_amplitude;
        Eigen::RowVectorXd ssa_factor;
        Eigen::RowVectorXd extinction_factor;
        BatchMatrix phase;
        BatchMatrix endpoint_source;

        void resize(int capacity) {
            ssa.resize(capacity);
            extinction.resize(capacity);
            unscaled_amplitude.resize(capacity);
            source_amplitude.resize(capacity);
            ssa_factor.resize(capacity);
            extinction_factor.resize(capacity);
            phase.resize(NSTOKES, capacity);
            endpoint_source.resize(NSTOKES, capacity);
        }
    };

    template <int NSTOKES> struct ExactIntegrationBlockScratch {
        using BatchMatrix =
            Eigen::Matrix<double, NSTOKES, Eigen::Dynamic, Eigen::RowMajor>;

        Eigen::RowVectorXd source_factor;
        Eigen::RowVectorXd source_factor_derivative;
        Eigen::RowVectorXd start_derivative_scale;
        Eigen::RowVectorXd end_derivative_scale;
        BatchMatrix integrated_value;
        BatchMatrix endpoint_quadrature;

        void resize(int capacity) {
            source_factor.resize(capacity);
            source_factor_derivative.resize(capacity);
            start_derivative_scale.resize(capacity);
            end_derivative_scale.resize(capacity);
            integrated_value.resize(NSTOKES, capacity);
            endpoint_quadrature.resize(NSTOKES, capacity);
        }
    };

    /** Batch equivalent of accumulate_exact_scattering_source. Values for a
     * fixed atmospheric/derivative coordinate are contiguous in wavelength so
     * Eigen can vectorize the common arithmetic. */
    template <int NSTOKES, int N>
    inline void accumulate_exact_scattering_source_block(
        const PhaseHandler<NSTOKES>& phase_handler, int threadidx, int losidx,
        int layeridx, const sasktran2::WavelengthBlock<N>& batch,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance,
        const Eigen::Ref<const Eigen::Matrix<double, 1, N, Eigen::RowMajor>>&
            solar_trans,
        const atmosphere::Atmosphere<NSTOKES>& atmosphere,
        SolarGeometryMatrix::InnerIterator solar_trans_iter,
        const Eigen::Ref<const Eigen::Matrix<double, 1, N, Eigen::RowMajor>>&
            derivative_scale,
        sasktran2::WavelengthBlockDual<NSTOKES>& target,
        ExactScatteringBlockScratch<NSTOKES>& scratch) {
        const auto& storage = atmosphere.storage();
        auto ssa = wavelength_head(scratch.ssa, batch);
        auto extinction = wavelength_head(scratch.extinction, batch);
        ssa.setZero();
        extinction.setZero();
        for (std::size_t index = 0; index < index_weights.size(); ++index) {
            const auto weight = index_weights[index];
            if (weight.second == 0.0) {
                continue;
            }
            for (int lane = 0; lane < batch.count; ++lane) {
                const int wavelength = batch.wavelength(lane);
                ssa(lane) +=
                    storage.ssa(weight.first, wavelength) * weight.second;
                extinction(lane) +=
                    storage.total_extinction(weight.first, wavelength) *
                    weight.second;
            }
        }

        auto unscaled_amplitude =
            wavelength_head(scratch.unscaled_amplitude, batch);
        auto source_amplitude =
            wavelength_head(scratch.source_amplitude, batch);
        unscaled_amplitude.array() = solar_trans.array() / (EIGEN_PI * 4);
        source_amplitude.array() =
            extinction.array() * ssa.array() * unscaled_amplitude.array();
        phase_handler.template scatter_and_accumulate_derivative_block<N>(
            threadidx, losidx, layeridx, index_weights, is_entrance, batch,
            source_amplitude, derivative_scale, target, scratch.phase);
        auto phase = wavelength_left_cols(scratch.phase, batch);
        auto endpoint_source =
            wavelength_left_cols(scratch.endpoint_source, batch);
        endpoint_source = phase;
        endpoint_source.array().rowwise() *= source_amplitude.array();

        if (target.derivative_size() > 0) {
            for (auto derivative = solar_trans_iter; derivative; ++derivative) {
                auto target_derivative =
                    target.derivative(derivative.index(), batch);
                target_derivative.array() -=
                    (endpoint_source.array().rowwise() *
                     derivative_scale.array()) *
                    derivative.value();
            }

            for (std::size_t index = 0; index < index_weights.size(); ++index) {
                const auto weight = index_weights[index];
                if (weight.second == 0.0) {
                    continue;
                }
                auto ssa_derivative = target.derivative(
                    atmosphere.ssa_deriv_start_index() + weight.first, batch);
                auto ssa_factor = wavelength_head(scratch.ssa_factor, batch);
                ssa_factor.array() = weight.second * derivative_scale.array() *
                                     extinction.array() *
                                     unscaled_amplitude.array();
                ssa_derivative.array() +=
                    phase.array().rowwise() * ssa_factor.array();

                auto extinction_derivative =
                    target.derivative(weight.first, batch);
                auto extinction_factor =
                    wavelength_head(scratch.extinction_factor, batch);
                extinction_factor.array() =
                    weight.second * derivative_scale.array() * ssa.array() *
                    unscaled_amplitude.array();
                extinction_derivative.array() +=
                    phase.array().rowwise() * extinction_factor.array();
            }
        }
    }

    template <typename S, int NSTOKES>
    class SingleScatterSource : public SourceTermInterface<NSTOKES> {
      private:
        static constexpr bool exact_transmission =
            std::is_same_v<S, SolarTransmissionExact>;
#ifdef SKTRAN_RUST_SUPPORT
        static constexpr bool compact_2d_table =
            std::is_same_v<S, SolarTransmissionTable2D>;
#else
        static constexpr bool compact_2d_table = false;
#endif
        static constexpr bool native_transmission_linearization =
            exact_transmission || compact_2d_table;

        std::shared_ptr<S> m_solar_transmission;
#ifdef SKTRAN_RUST_SUPPORT
        std::shared_ptr<SolarTransmissionTable2D> m_shared_solar_table_2d;
        const sasktran2::raytracing::RustRayTracer2D* m_raytracer_2d = nullptr;
#endif
        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere =
            nullptr;
        bool m_native_volume_linearization_active = true;
        std::uint64_t m_atmosphere_instance_id = 0;
        std::uint64_t m_atmosphere_volume_revision = 0;
        bool m_has_atmosphere_volume_revision = false;

        Eigen::MatrixXd m_geometry_matrix;
        SolarGeometryMatrix m_geometry_sparse;
        std::vector<bool> m_ground_hit_flag;

        std::vector<Eigen::VectorXd> m_solar_trans;
        using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic,
                                             Eigen::Dynamic, Eigen::RowMajor>;
        std::vector<RowMajorMatrix> m_solar_trans_batch;
        SolarTableInterpolation m_solar_interpolation;
        std::vector<Eigen::Vector3d> m_solar_propagation_directions;
        std::vector<Eigen::VectorXd> m_solar_table_product;
        std::vector<Eigen::VectorXd> m_solar_trans_jvp;
        mutable std::vector<Eigen::VectorXd> m_solar_endpoint_cotangent;
        mutable Eigen::VectorXd m_solar_endpoint_cotangent_sum;
        mutable Eigen::VectorXd m_solar_table_cotangent;
        int m_wavelength_batch_capacity = 1;
        std::vector<int> m_active_wavelength_block_start;
        std::vector<int> m_active_wavelength_block_count;
        std::vector<std::vector<int>> m_index_map;

        PhaseHandler<NSTOKES> m_phase_handler;

        // [los][layer][solar endpoint][geometry endpoint], where endpoint 0 is
        // the layer exit and endpoint 1 is the layer entrance.
        using ActiveDerivativeIndices =
            std::array<std::array<std::vector<int>, 2>, 2>;
        std::vector<std::vector<ActiveDerivativeIndices>>
            m_active_derivative_indices;
        const std::vector<sasktran2::raytracing::TracedRay>* m_traced_rays =
            nullptr;

        mutable std::vector<
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>>
            m_start_source_cache;
        mutable std::vector<
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>>
            m_end_source_cache;
        mutable std::vector<std::array<ExactScatteringBlockScratch<NSTOKES>, 2>>
            m_batch_source_cache;
        mutable std::vector<ExactIntegrationBlockScratch<NSTOKES>>
            m_batch_integration_cache;

        const Geometry& m_geometry;
        const Geometry1D* m_geometry_1d = nullptr;
        const Geometry2D* m_geometry_2d = nullptr;
        const sasktran2::Config* m_config;

        std::vector<bool> m_los_ground_is_hit;
        std::vector<sasktran2::raytracing::LayerGeometry> m_los_end_layers;
        std::vector<std::vector<std::pair<int, double>>>
            m_los_surface_interpolation_weights;

        void initialize_active_derivative_indices();
        void initialize_atmosphere_impl(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
            bool materialized_derivative_storage);

        void initialize_fixed_dispatch() {
            if constexpr (std::is_same_v<S, SolarTransmissionExact>) {
                this->template set_fixed_integrated_source_dispatch<1>(
                    &fixed_integrated_source_dispatch<1>);
                this->template set_fixed_integrated_source_dispatch<4>(
                    &fixed_integrated_source_dispatch<4>);
            }
        }

        void integrated_source_constant(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::LayerGeometry& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            sasktran2::WavelengthBlockLaneDualView<NSTOKES, 1>& source,
            typename SourceTermInterface<NSTOKES>::IntegrationDirection
                direction) const;

        void endpoint_source_jvp(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int solar_index,
            const sasktran2::raytracing::GridWeightStencilView& weights,
            bool is_entrance, Eigen::Ref<const Eigen::VectorXd> native_tangent,
            sasktran2::RadianceJVP<NSTOKES>& result) const;

        Eigen::Vector<double, NSTOKES> endpoint_source_vjp(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, int solar_index,
            const sasktran2::raytracing::GridWeightStencilView& weights,
            bool is_entrance, const Eigen::Vector<double, NSTOKES>& cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const;

        double solar_transmission_value(int wavelidx, int threadidx,
                                        int solar_index) const;
        double solar_transmission_tangent(
            int wavelidx, int threadidx, int solar_index,
            Eigen::Ref<const Eigen::VectorXd> native_tangent) const;
        bool ground_scattering_geometry(int losidx, double& mu_in,
                                        double& mu_out, double& phi_diff) const;

#ifdef SKTRAN_RUST_SUPPORT
        static std::shared_ptr<S> make_2d_transmission(
            const Geometry2D& geometry,
            const sasktran2::raytracing::RustRayTracer2D& raytracer,
            const std::shared_ptr<SolarTransmissionTable2D>& shared_table) {
            if constexpr (compact_2d_table) {
                return shared_table != nullptr
                           ? shared_table
                           : std::make_shared<SolarTransmissionTable2D>(
                                 geometry, raytracer);
            } else {
                if constexpr (exact_transmission) {
                    return std::make_shared<S>(geometry, raytracer);
                } else {
                    return {};
                }
            }
        }
#endif

      public:
        template <
            typename T = S,
            std::enable_if_t<std::is_same_v<T, SolarTransmissionExact> ||
                                 std::is_same_v<T, SolarTransmissionTable>,
                             int> = 0>
        SingleScatterSource(
            const Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer)
            : m_solar_transmission(std::make_shared<S>(geometry, raytracer)),
              m_geometry(geometry), m_geometry_1d(&geometry),
              m_phase_handler(geometry) {
            initialize_fixed_dispatch();
        };

#ifdef SKTRAN_RUST_SUPPORT
        template <
            typename T = S,
            std::enable_if_t<std::is_same_v<T, SolarTransmissionExact> ||
                                 std::is_same_v<T, SolarTransmissionTable2D>,
                             int> = 0>
        SingleScatterSource(
            const Geometry2D& geometry,
            const sasktran2::raytracing::RustRayTracer2D& raytracer,
            std::shared_ptr<SolarTransmissionTable2D> shared_table = nullptr)
            : m_solar_transmission(
                  make_2d_transmission(geometry, raytracer, shared_table)),
              m_shared_solar_table_2d(std::move(shared_table)),
              m_raytracer_2d(&raytracer), m_geometry(geometry),
              m_geometry_2d(&geometry), m_phase_handler(geometry) {
            if constexpr (compact_2d_table) {
                m_shared_solar_table_2d = m_solar_transmission;
            }
            initialize_fixed_dispatch();
        };
#endif

        void initialize_config(const sasktran2::Config& config) override;

        /** Here the single scatter source term initializes the internal solar
         * transmission object, usually this involves tracing the required rays
         * and setting up any internal matrices.
         *
         *  @param internal_viewing Information on the internal viewing
         * geometry, los_rays and flux observers
         */
        void initialize_geometry(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing) override;

        /**
         *
         */
        void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere)
            override;
        void initialize_atmosphere_native(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere)
            override;

        void set_wavelength_block_capacity(int block_capacity) override {
            if (block_capacity < 1) {
                throw std::invalid_argument(
                    "Single scatter wavelength block capacity must be "
                    "positive");
            }
            if (block_capacity != m_wavelength_batch_capacity) {
                std::fill(m_active_wavelength_block_count.begin(),
                          m_active_wavelength_block_count.end(), 0);
                m_has_atmosphere_volume_revision = false;
            }
            m_wavelength_batch_capacity = block_capacity;
        }

        /** Triggers an internal calculation of the source term.  This method is
         * called at the beginning of each 'wavelength' calculation.
         *
         * @param wavelidx Index of the wavelength being calculated
         */
        int maximum_wavelength_block_size() const override {
            return std::is_same_v<S, SolarTransmissionExact>
                       ? std::numeric_limits<int>::max()
                       : 1;
        }

      private:
        void calculate_single(int wavelidx, int threadidx);

        void initialize_wavelength_blocks(int block_size);

        template <int N>
        void calculate_block(const sasktran2::WavelengthBlock<N>& block,
                             int threadidx);

      public:
        void calculate(const sasktran2::WavelengthBlock<>& block,
                       int threadidx) override {
            if constexpr (exact_transmission) {
                sasktran2::dispatch_wavelength_block(
                    block, [&](const auto& fixed_block) {
                        if constexpr (std::decay_t<
                                          decltype(fixed_block)>::static_size ==
                                      1) {
                            if (m_wavelength_batch_capacity == 1) {
                                calculate_single(fixed_block.start, threadidx);
                            } else {
                                calculate_block(fixed_block, threadidx);
                            }
                        } else {
                            calculate_block(fixed_block, threadidx);
                        }
                    });
            } else {
                calculate_single(block.start, threadidx);
            }
        }

        void calculate_jvp(
            const sasktran2::WavelengthBlock<>& block, int threadidx,
            Eigen::Ref<const Eigen::VectorXd> native_tangent) override;
        void calculate_vjp(const sasktran2::WavelengthBlock<>& block,
                           int threadidx) override;
        void finalize_vjp(
            const sasktran2::WavelengthBlock<>& block, int threadidx,
            Eigen::Ref<Eigen::MatrixXd> native_gradient) const override;

        /** Calculates the integrated source term for a given layer.
         *
         * @param losidx Raw index pointing to the ray that was previously
         * passed in initialize_geometry
         * @param layeridx Raw index pointing to the layer that was previosuly
         * passed in initialize_geometry
         * @param layer The layer that we are integrating over
         * @param source The returned source term
         */
      private:
        void integrated_source_single(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            sasktran2::WavelengthBlockLaneDualView<NSTOKES, 1>& source,
            typename SourceTermInterface<
                NSTOKES>::IntegrationDirection direction =
                SourceTermInterface<NSTOKES>::IntegrationDirection::none) const;

        template <int N>
        void integrated_source_block(
            const sasktran2::WavelengthBlock<N>& batch, int losidx,
            int layeridx, int wavel_threadidx, int threadidx,
            const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            sasktran2::WavelengthBlockDual<NSTOKES>& source,
            typename SourceTermInterface<
                NSTOKES>::IntegrationDirection direction =
                SourceTermInterface<NSTOKES>::IntegrationDirection::none) const;

        template <int N>
        void integrated_source_typed(
            const sasktran2::WavelengthBlock<N>& block, int losidx,
            int layeridx, int wavel_threadidx, int threadidx,
            const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            sasktran2::WavelengthBlockDual<NSTOKES>& source,
            typename SourceTermInterface<NSTOKES>::IntegrationDirection
                direction) const {
            if constexpr (N == 1) {
                if (m_wavelength_batch_capacity == 1) {
                    sasktran2::WavelengthBlockLaneDualView<NSTOKES, 1>
                        source_lane(source, 0);
                    integrated_source_single(block.start, losidx, layeridx,
                                             wavel_threadidx, threadidx, layer,
                                             entrance_weights, exit_weights,
                                             shell_od, source_lane, direction);
                    return;
                }
            }
            integrated_source_block(block, losidx, layeridx, wavel_threadidx,
                                    threadidx, layer, entrance_weights,
                                    exit_weights, shell_od, source, direction);
        }

        template <int N>
        static void fixed_integrated_source_dispatch(
            const SourceTermInterface<NSTOKES>& source_term,
            const sasktran2::WavelengthBlock<N>& block, int losidx,
            int layeridx, int wavel_threadidx, int threadidx,
            const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            sasktran2::WavelengthBlockDual<NSTOKES>& source,
            typename SourceTermInterface<NSTOKES>::IntegrationDirection
                direction) {
            static_cast<const SingleScatterSource&>(source_term)
                .integrated_source_typed(block, losidx, layeridx,
                                         wavel_threadidx, threadidx, layer,
                                         entrance_weights, exit_weights,
                                         shell_od, source, direction);
        }

      public:
        void integrated_source(
            const sasktran2::WavelengthBlock<>& block, int losidx, int layeridx,
            int wavel_threadidx, int threadidx,
            const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            sasktran2::WavelengthBlockDual<NSTOKES>& source,
            typename SourceTermInterface<NSTOKES>::IntegrationDirection
                direction =
                    SourceTermInterface<NSTOKES>::IntegrationDirection::none)
            const override {
            if constexpr (std::is_same_v<S, SolarTransmissionExact>) {
                sasktran2::dispatch_wavelength_block(
                    block, [&](const auto& fixed_block) {
                        integrated_source_typed(
                            fixed_block, losidx, layeridx, wavel_threadidx,
                            threadidx, layer, entrance_weights, exit_weights,
                            shell_od, source, direction);
                    });
            } else {
                sasktran2::WavelengthBlockLaneDualView<NSTOKES, 1> source_lane(
                    source, 0);
                integrated_source_single(block.start, losidx, layeridx,
                                         wavel_threadidx, threadidx, layer,
                                         entrance_weights, exit_weights,
                                         shell_od, source_lane, direction);
            }
        }

        bool supports_geometry_dimension(int dimension) const override {
            return dimension == 1 ||
                   (dimension == 2 && m_geometry_2d != nullptr);
        }

        bool supports_sparse_derivative_tracking() const override {
            return exact_transmission;
        }

        bool supports_linearization(
            sasktran2::LinearizationMode mode) const override {
            if constexpr (exact_transmission) {
                return true;
            }
            if constexpr (compact_2d_table) {
                return mode != sasktran2::LinearizationMode::Jacobian;
            }
            return mode == sasktran2::LinearizationMode::Jacobian;
        }

        void end_of_ray_source_jvp(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            Eigen::Ref<const Eigen::VectorXd> native_tangent,
            sasktran2::RadianceJVP<NSTOKES>& source) const override;

        void integrated_source_jvp(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            Eigen::Ref<const Eigen::VectorXd> native_tangent,
            sasktran2::RadianceJVP<NSTOKES>& source) const override;

        /** Native exact single scattering is additive at the end of the ray
         * and within each layer, so these implementations leave the mutable
         * cotangent unchanged. */
        void end_of_ray_source_vjp(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            const Eigen::Vector<double, NSTOKES>& value_before,
            Eigen::Vector<double, NSTOKES>& cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const override;

        void integrated_source_vjp(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            const Eigen::Vector<double, NSTOKES>& value_before,
            Eigen::Vector<double, NSTOKES>& cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const override;

        void append_end_of_ray_active_derivatives(
            int losidx, std::vector<int>& derivative_indices) const override;

        void append_interior_active_derivatives(
            int losidx, int layeridx,
            std::vector<int>& derivative_indices) const override;

        /** Calculates the source term at the end of the ray.  Common examples
         * of this are ground scattering, ground emission, or the solar radiance
         * if looking directly at the sun.
         *
         * @param losidx Raw index pointing to the ray that was previously
         * passed in initialize_geometry
         * @param surface The surface object
         * @param source The returned source term
         */
      private:
        void end_of_ray_source_single(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            sasktran2::WavelengthBlockLaneDualView<NSTOKES, 1>& source) const;

        template <int N>
        void end_of_ray_source_block(
            const sasktran2::WavelengthBlock<N>& batch, int losidx,
            int wavel_threadidx, int threadidx,
            sasktran2::WavelengthBlockDual<NSTOKES>& source) const;

      public:
        void end_of_ray_source(
            const sasktran2::WavelengthBlock<>& block, int losidx,
            int wavel_threadidx, int threadidx,
            sasktran2::WavelengthBlockDual<NSTOKES>& source) const override {
            if constexpr (std::is_same_v<S, SolarTransmissionExact>) {
                sasktran2::dispatch_wavelength_block(
                    block, [&](const auto& fixed_block) {
                        if constexpr (std::decay_t<
                                          decltype(fixed_block)>::static_size ==
                                      1) {
                            if (m_wavelength_batch_capacity == 1) {
                                sasktran2::WavelengthBlockLaneDualView<NSTOKES,
                                                                       1>
                                    source_lane(source, 0);
                                end_of_ray_source_single(
                                    fixed_block.start, losidx, wavel_threadidx,
                                    threadidx, source_lane);
                            } else {
                                end_of_ray_source_block(fixed_block, losidx,
                                                        wavel_threadidx,
                                                        threadidx, source);
                            }
                        } else {
                            end_of_ray_source_block(fixed_block, losidx,
                                                    wavel_threadidx, threadidx,
                                                    source);
                        }
                    });
            } else {
                sasktran2::WavelengthBlockLaneDualView<NSTOKES, 1> source_lane(
                    source, 0);
                end_of_ray_source_single(block.start, losidx, wavel_threadidx,
                                         threadidx, source_lane);
            }
        }

        /**
         * @brief Not used for the Single Scatter source.
         *
         * @param wavelidx
         * @param losidx
         * @param wavel_threadidx
         * @param threadidx
         * @param source
         */
        void start_of_ray_source(
            const sasktran2::WavelengthBlock<>&, int, int, int,
            sasktran2::WavelengthBlockDual<NSTOKES>&) const override {}

        void start_of_ray_source_jvp(
            int, int, int, int, Eigen::Ref<const Eigen::VectorXd>,
            sasktran2::RadianceJVP<NSTOKES>&) const override {}

        void
        start_of_ray_source_vjp(int, int, int, int,
                                const Eigen::Vector<double, NSTOKES>&,
                                Eigen::Vector<double, NSTOKES>&,
                                Eigen::Ref<Eigen::VectorXd>) const override {}
    };

    template <int NSTOKES>
    class OccultationSource : public SourceTermInterface<NSTOKES> {
      private:
        std::vector<bool> m_ground_is_hit;

      public:
        void initialize_config(const sasktran2::Config& config) override;

        /** Initializes any geometry information that is required for
         * calculating the source term.  This method is called after the line of
         * sight rays ar traced.
         *
         * @param internal_viewing Information on the internal viewing geometry,
         * los_rays and flux observers
         */
        void initialize_geometry(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing) override;

        /**
         *
         */
        void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere)
            override;

        /** Triggers an internal calculation of the source term.  This method is
         * called at the beginning of each 'wavelength' calculation.
         *
         * @param wavelidx Index of the wavelength being calculated
         */
        int maximum_wavelength_block_size() const override {
            return std::numeric_limits<int>::max();
        }

        void calculate(const sasktran2::WavelengthBlock<>&, int) override {}

        /** Calculates the integrated source term for a given layer.
         *
         * @param losidx Raw index pointing to the ray that was previously
         * passed in initialize_geometry
         * @param layeridx Raw index pointing to the layer that was previosuly
         * passed in initialize_geometry
         * @param layer The layer that we are integrating over
         * @param source The returned source term
         */
        void integrated_source(
            const sasktran2::WavelengthBlock<>&, int, int, int, int,
            const sasktran2::raytracing::TracedLayer&,
            const sasktran2::raytracing::GridWeightStencilView&,
            const sasktran2::raytracing::GridWeightStencilView&,
            const sasktran2::WavelengthBlockODView&,
            sasktran2::WavelengthBlockDual<NSTOKES>&,
            typename SourceTermInterface<NSTOKES>::IntegrationDirection =
                SourceTermInterface<NSTOKES>::IntegrationDirection::none)
            const override {}

        /** Calculates the source term at the end of the ray.  Common examples
         * of this are ground scattering, ground emission, or the solar radiance
         * if looking directly at the sun.
         *
         * @param wavelidx Raw index for the wavelength we are calculating
         * @param losidx Raw index pointing to the ray that was previously
         * passed in initialize_geometry
         * @param source The returned source term
         */
        void end_of_ray_source(
            const sasktran2::WavelengthBlock<>& block, int losidx,
            int wavel_threadidx, int threadidx,
            sasktran2::WavelengthBlockDual<NSTOKES>& source) const override;

        /**
         * @brief Not used for the occultation source.
         *
         * @param wavelidx
         * @param losidx
         * @param wavel_threadidx
         * @param threadidx
         * @param source
         */
        void start_of_ray_source(
            const sasktran2::WavelengthBlock<>&, int, int, int,
            sasktran2::WavelengthBlockDual<NSTOKES>&) const override {}

        bool has_interior_source() const override { return false; }
    };

} // namespace sasktran2::solartransmission
