#include "sasktran2/geometry.h"
#include <sasktran2/solartransmission.h>
#include <type_traits>
#include <unordered_map>

namespace sasktran2::solartransmission {

    void SolarTableInterpolation::clear() {
        std::vector<std::uint32_t>().swap(m_outer);
        std::vector<std::uint32_t>().swap(m_inner);
        std::vector<std::uint8_t>().swap(m_row_counts);
        std::vector<std::uint16_t>().swap(m_relative_inner);
        std::vector<std::uint32_t>().swap(m_row_bases);
        std::vector<std::uint16_t>().swap(m_row_patterns16);
        std::vector<std::uint32_t>().swap(m_row_patterns32);
        std::vector<std::uint32_t>().swap(m_column_patterns);
        std::vector<double>().swap(m_values);
        m_rows = 0;
        m_cols = 0;
        m_next_row = 0;
        m_compact_rows = false;
        m_relative_indices = false;
        m_interned_indices = false;
        m_finalized = false;
    }

    void SolarTableInterpolation::initialize(Eigen::Index rows,
                                             Eigen::Index cols,
                                             Eigen::Index capacity) {
        if (rows < 0 || cols < 0 || capacity < 0 ||
            cols > static_cast<Eigen::Index>(
                       std::numeric_limits<std::uint32_t>::max())) {
            throw std::invalid_argument(
                "Invalid compact solar interpolation dimensions");
        }
        m_rows = rows;
        m_cols = cols;
        m_next_row = 0;
        m_compact_rows = false;
        m_relative_indices = false;
        m_interned_indices = false;
        m_finalized = false;
        std::vector<std::uint8_t>().swap(m_row_counts);
        std::vector<std::uint16_t>().swap(m_relative_inner);
        std::vector<std::uint32_t>().swap(m_row_bases);
        std::vector<std::uint16_t>().swap(m_row_patterns16);
        std::vector<std::uint32_t>().swap(m_row_patterns32);
        std::vector<std::uint32_t>().swap(m_column_patterns);
        m_outer.assign(static_cast<std::size_t>(rows) + 1, 0);
        m_inner.clear();
        m_values.clear();
        m_inner.reserve(static_cast<std::size_t>(capacity));
        m_values.reserve(static_cast<std::size_t>(capacity));
    }

    void SolarTableInterpolation::append_row(
        const std::vector<std::pair<int, double>>& interpolation_weights) {
        if (m_next_row >= m_rows) {
            throw std::logic_error(
                "Too many rows appended to compact solar interpolation");
        }
        m_outer[static_cast<std::size_t>(m_next_row)] =
            static_cast<std::uint32_t>(m_values.size());
        for (const auto& [column, value] : interpolation_weights) {
            if (column < 0 || column >= m_cols || !std::isfinite(value)) {
                throw std::invalid_argument(
                    "Invalid compact solar interpolation entry");
            }
            if (value == 0.0) {
                continue;
            }
            if (m_values.size() == std::numeric_limits<std::uint32_t>::max()) {
                throw std::length_error(
                    "Compact solar interpolation has too many entries");
            }
            m_inner.push_back(static_cast<std::uint32_t>(column));
            m_values.push_back(value);
        }
        ++m_next_row;
    }

    bool SolarTableInterpolation::try_intern_column_patterns(
        std::size_t maximum_index_bytes) {
        const auto rows = static_cast<std::size_t>(m_rows);
        if (rows * sizeof(std::uint16_t) >= maximum_index_bytes) {
            return false;
        }
        struct Pattern {
            std::uint32_t begin;
            std::uint32_t count;
            std::uint32_t id;
        };
        std::unordered_multimap<std::uint64_t, Pattern> patterns;
        patterns.reserve(std::min<std::size_t>(rows, 65536));
        std::vector<std::uint16_t> row_ids16;
        std::vector<std::uint32_t> row_ids32;
        std::vector<std::uint32_t> descriptors;
        std::vector<std::uint32_t> columns;
        row_ids16.reserve(rows);
        descriptors.reserve(std::min<std::size_t>(rows, 65536));
        columns.reserve(std::min<std::size_t>(m_inner.size(), 65536));
        bool wide_ids = false;
        for (Eigen::Index row = 0; row < m_rows; ++row) {
            const auto begin = m_outer[static_cast<std::size_t>(row)];
            const auto end = m_outer[static_cast<std::size_t>(row) + 1];
            const auto count = end - begin;
            std::uint64_t hash = 14695981039346656037ULL;
            hash = (hash ^ count) * 1099511628211ULL;
            for (std::uint32_t entry = begin; entry < end; ++entry) {
                hash = (hash ^ m_inner[entry]) * 1099511628211ULL;
            }
            bool found = false;
            std::uint32_t id = 0;
            const auto matching_hashes = patterns.equal_range(hash);
            for (auto pattern = matching_hashes.first;
                 pattern != matching_hashes.second; ++pattern) {
                const auto& representative = pattern->second;
                if (representative.count != count) {
                    continue;
                }
                bool equal = true;
                for (std::uint32_t entry = 0; entry < count; ++entry) {
                    if (m_inner[begin + entry] !=
                        m_inner[representative.begin + entry]) {
                        equal = false;
                        break;
                    }
                }
                if (equal) {
                    found = true;
                    id = representative.id;
                    break;
                }
            }
            if (!found) {
                // The descriptor's low four bits store the count; its high
                // 28 bits address the shared ordered column array.
                constexpr std::size_t maximum_dictionary_entries =
                    std::size_t{1} << 28;
                if (count > 15 ||
                    columns.size() >= maximum_dictionary_entries ||
                    columns.size() + count > maximum_dictionary_entries) {
                    return false;
                }
                id = static_cast<std::uint32_t>(descriptors.size());
                const bool need_wide_ids =
                    id > std::numeric_limits<std::uint16_t>::max();
                const auto id_bytes = need_wide_ids ? sizeof(std::uint32_t)
                                                    : sizeof(std::uint16_t);
                const auto retained_bytes =
                    rows * id_bytes +
                    (descriptors.size() + 1 + columns.size() + count) *
                        sizeof(std::uint32_t);
                // Stop growing the temporary dictionary as soon as even its
                // current unique payload cannot improve the existing format.
                if (retained_bytes >= maximum_index_bytes) {
                    return false;
                }
                if (need_wide_ids && !wide_ids) {
                    row_ids32.reserve(rows);
                    row_ids32.insert(row_ids32.end(), row_ids16.begin(),
                                     row_ids16.end());
                    std::vector<std::uint16_t>().swap(row_ids16);
                    wide_ids = true;
                }
                descriptors.push_back(
                    (static_cast<std::uint32_t>(columns.size()) << 4) | count);
                for (std::uint32_t entry = begin; entry < end; ++entry) {
                    columns.push_back(m_inner[entry]);
                }
                patterns.emplace(hash, Pattern{begin, count, id});
            }
            if (wide_ids) {
                row_ids32.push_back(id);
            } else {
                row_ids16.push_back(static_cast<std::uint16_t>(id));
            }
        }
        row_ids16.shrink_to_fit();
        row_ids32.shrink_to_fit();
        descriptors.shrink_to_fit();
        columns.shrink_to_fit();
        const std::size_t allocated_bytes =
            row_ids16.capacity() * sizeof(std::uint16_t) +
            row_ids32.capacity() * sizeof(std::uint32_t) +
            (descriptors.capacity() + columns.capacity()) *
                sizeof(std::uint32_t);
        if (allocated_bytes >= maximum_index_bytes) {
            return false;
        }
        // Only complete arrays are committed. Original row weights and their
        // floating-point operation sequence remain separate and unchanged.
        m_row_patterns16.swap(row_ids16);
        m_row_patterns32.swap(row_ids32);
        m_column_patterns.swap(descriptors);
        m_inner.swap(columns);
        std::vector<std::uint32_t>().swap(m_outer);
        return true;
    }

    void SolarTableInterpolation::finalize() {
        if (m_finalized) {
            return;
        }
        if (m_next_row != m_rows) {
            throw std::logic_error(
                "Compact solar interpolation has missing rows");
        }
        m_outer[static_cast<std::size_t>(m_rows)] =
            static_cast<std::uint32_t>(m_values.size());
        m_inner.shrink_to_fit();
        m_values.shrink_to_fit();
        const auto original_storage_bytes = storage_bytes();
        std::uint32_t maximum_row_span = 0;
        std::uint32_t maximum_row_nonzeros = 0;
        std::size_t wide_span_rows = 0;
        for (Eigen::Index row = 0; row < m_rows; ++row) {
            const auto begin = m_outer[static_cast<std::size_t>(row)];
            const auto end = m_outer[static_cast<std::size_t>(row) + 1];
            maximum_row_nonzeros = std::max(maximum_row_nonzeros, end - begin);
            if (begin == end) {
                continue;
            }
            std::uint32_t minimum = m_inner[begin];
            std::uint32_t maximum = minimum;
            for (std::uint32_t entry = begin + 1; entry < end; ++entry) {
                minimum = std::min(minimum, m_inner[entry]);
                maximum = std::max(maximum, m_inner[entry]);
            }
            const auto span = maximum - minimum;
            maximum_row_span = std::max(maximum_row_span, span);
            wide_span_rows += span > std::numeric_limits<std::uint16_t>::max();
        }
        // Every relative row must fit. Keep the wide format when the added
        // bases would outweigh the saved column-index bytes.
        m_relative_indices =
            wide_span_rows == 0 &&
            m_inner.size() * sizeof(std::uint16_t) >
                static_cast<std::size_t>(m_rows) * sizeof(std::uint32_t);
        const std::size_t best_existing_index_bytes =
            (maximum_row_nonzeros <= std::numeric_limits<std::uint8_t>::max()
                 ? static_cast<std::size_t>(m_rows) * sizeof(std::uint8_t)
                 : (static_cast<std::size_t>(m_rows) + 1) *
                       sizeof(std::uint32_t)) +
            (m_relative_indices
                 ? m_inner.size() * sizeof(std::uint16_t) +
                       static_cast<std::size_t>(m_rows) * sizeof(std::uint32_t)
                 : m_inner.size() * sizeof(std::uint32_t));
        m_interned_indices =
            maximum_row_nonzeros <= 15 &&
            try_intern_column_patterns(best_existing_index_bytes);
        if (m_interned_indices) {
            m_relative_indices = false;
        }
        if (m_relative_indices) {
            m_row_bases.resize(static_cast<std::size_t>(m_rows));
            m_relative_inner.resize(m_inner.size());
            for (Eigen::Index row = 0; row < m_rows; ++row) {
                const auto begin = m_outer[static_cast<std::size_t>(row)];
                const auto end = m_outer[static_cast<std::size_t>(row) + 1];
                const auto base =
                    begin == end ? 0
                                 : *std::min_element(m_inner.begin() + begin,
                                                     m_inner.begin() + end);
                m_row_bases[static_cast<std::size_t>(row)] = base;
                for (std::uint32_t entry = begin; entry < end; ++entry) {
                    m_relative_inner[entry] =
                        static_cast<std::uint16_t>(m_inner[entry] - base);
                }
            }
            std::vector<std::uint32_t>().swap(m_inner);
        }
        // Products consume complete rows in their original sequence; no
        // random row lookup needs the prefix array after construction.
        m_compact_rows =
            !m_interned_indices &&
            maximum_row_nonzeros <= std::numeric_limits<std::uint8_t>::max();
        if (m_compact_rows) {
            m_row_counts.resize(static_cast<std::size_t>(m_rows));
            for (Eigen::Index row = 0; row < m_rows; ++row) {
                m_row_counts[static_cast<std::size_t>(row)] =
                    static_cast<std::uint8_t>(
                        m_outer[static_cast<std::size_t>(row) + 1] -
                        m_outer[static_cast<std::size_t>(row)]);
            }
            std::vector<std::uint32_t>().swap(m_outer);
        }
        m_finalized = true;
        if (std::getenv("SASKTRAN2_PROFILE_MEMORY") != nullptr) {
            std::fprintf(
                stderr,
                "SASKTRAN2_MEMORY "
                "{\"kind\":\"solar_table_interpolation_memory\","
                "\"owner\":\"%p\",\"rows\":%lld,\"columns\":%lld,"
                "\"nonzeros\":%zu,\"maximum_row_span\":%u,"
                "\"wide_span_rows\":%zu,\"maximum_row_nonzeros\":%u,"
                "\"relative_column_indices\":%s,\"compact_row_counts\":%s,"
                "\"interned_column_patterns\":%s,\"pattern_id_bits\":%d,"
                "\"unique_patterns\":%zu,\"dictionary_descriptor_bytes\":%zu,"
                "\"row_pattern_id_bytes\":%zu,\"dictionary_column_bytes\":%zu,"
                "\"best_existing_index_bytes\":%zu,"
                "\"retained_index_bytes\":%zu,"
                "\"interned_index_bytes_saved\":%zu,"
                "\"outer_size\":%zu,\"outer_capacity\":%zu,"
                "\"inner_size\":%zu,\"inner_capacity\":%zu,"
                "\"relative_inner_size\":%zu,\"relative_inner_capacity\":%zu,"
                "\"row_bases_size\":%zu,\"row_bases_capacity\":%zu,"
                "\"row_counts_size\":%zu,\"row_counts_capacity\":%zu,"
                "\"value_size\":%zu,\"value_capacity\":%zu,"
                "\"original_allocated_bytes\":%zu,\"allocated_bytes\":%zu,"
                "\"column_index_bytes_saved\":%zu,\"row_offset_bytes_saved\":%"
                "zu}\n",
                static_cast<const void*>(this), static_cast<long long>(m_rows),
                static_cast<long long>(m_cols), m_values.size(),
                maximum_row_span, wide_span_rows, maximum_row_nonzeros,
                m_relative_indices ? "true" : "false",
                m_compact_rows ? "true" : "false",
                m_interned_indices ? "true" : "false", pattern_id_bits(),
                m_column_patterns.size(),
                m_column_patterns.capacity() * sizeof(std::uint32_t),
                m_row_patterns16.capacity() * sizeof(std::uint16_t) +
                    m_row_patterns32.capacity() * sizeof(std::uint32_t),
                m_interned_indices ? m_inner.capacity() * sizeof(std::uint32_t)
                                   : 0,
                best_existing_index_bytes,
                storage_bytes() - m_values.capacity() * sizeof(double),
                m_interned_indices ? best_existing_index_bytes -
                                         (storage_bytes() -
                                          m_values.capacity() * sizeof(double))
                                   : 0,
                m_outer.size(), m_outer.capacity(), m_inner.size(),
                m_inner.capacity(), m_relative_inner.size(),
                m_relative_inner.capacity(), m_row_bases.size(),
                m_row_bases.capacity(), m_row_counts.size(),
                m_row_counts.capacity(), m_values.size(), m_values.capacity(),
                original_storage_bytes, storage_bytes(),
                m_relative_indices
                    ? m_values.size() * sizeof(std::uint16_t) -
                          m_row_bases.size() * sizeof(std::uint32_t)
                    : 0,
                m_compact_rows ? (static_cast<std::size_t>(m_rows) + 1) *
                                         sizeof(std::uint32_t) -
                                     m_row_counts.size() * sizeof(std::uint8_t)
                               : 0);
        }
    }

    void SolarTableInterpolation::apply(
        Eigen::Ref<const Eigen::VectorXd> table_values,
        Eigen::Ref<Eigen::VectorXd> endpoint_values) const {
        if (table_values.size() != m_cols || endpoint_values.size() != m_rows) {
            throw std::invalid_argument(
                "Invalid compact solar interpolation product dimensions");
        }
        if (m_interned_indices) {
            const auto apply_patterns = [&](const auto& row_ids) {
                std::uint32_t position = 0;
                const double* values = table_values.data();
                for (Eigen::Index row = 0; row < m_rows; ++row) {
                    const auto descriptor =
                        m_column_patterns[row_ids[static_cast<std::size_t>(
                            row)]];
                    std::uint32_t column_position = descriptor >> 4;
                    const std::uint32_t end = position + (descriptor & 15);
                    double result = 0.0;
                    for (; position < end; ++position, ++column_position) {
                        result += m_values[position] *
                                  values[m_inner[column_position]];
                    }
                    endpoint_values(row) = result;
                }
            };
            if (m_row_patterns32.empty()) {
                apply_patterns(m_row_patterns16);
            } else {
                apply_patterns(m_row_patterns32);
            }
            return;
        }
        const auto apply_rows = [&](const auto& indices, const auto& rows,
                                    auto relative) {
            constexpr bool compact_rows = std::is_same_v<
                typename std::decay_t<decltype(rows)>::value_type,
                std::uint8_t>;
            std::uint32_t position = 0;
            for (Eigen::Index row = 0; row < m_rows; ++row) {
                std::uint32_t begin;
                std::uint32_t end;
                if constexpr (compact_rows) {
                    begin = position;
                    end = begin + rows[static_cast<std::size_t>(row)];
                    position = end;
                } else {
                    begin = rows[static_cast<std::size_t>(row)];
                    end = rows[static_cast<std::size_t>(row) + 1];
                }
                const double* row_values = table_values.data();
                if constexpr (decltype(relative)::value) {
                    if (begin != end) {
                        row_values +=
                            m_row_bases[static_cast<std::size_t>(row)];
                    }
                }
                double result = 0.0;
                for (std::uint32_t entry = begin; entry < end; ++entry) {
                    result += m_values[entry] * row_values[indices[entry]];
                }
                endpoint_values(row) = result;
            }
        };
        if (m_relative_indices) {
            if (m_compact_rows) {
                apply_rows(m_relative_inner, m_row_counts, std::true_type{});
            } else {
                apply_rows(m_relative_inner, m_outer, std::true_type{});
            }
        } else if (m_compact_rows) {
            apply_rows(m_inner, m_row_counts, std::false_type{});
        } else {
            apply_rows(m_inner, m_outer, std::false_type{});
        }
    }

    void SolarTableInterpolation::apply_transpose(
        Eigen::Ref<const Eigen::VectorXd> endpoint_values,
        Eigen::Ref<Eigen::VectorXd> table_values) const {
        if (endpoint_values.size() != m_rows || table_values.size() != m_cols) {
            throw std::invalid_argument(
                "Invalid compact solar interpolation transpose dimensions");
        }
        table_values.setZero();
        if (m_interned_indices) {
            const auto apply_patterns = [&](const auto& row_ids) {
                std::uint32_t position = 0;
                double* values = table_values.data();
                for (Eigen::Index row = 0; row < m_rows; ++row) {
                    const auto descriptor =
                        m_column_patterns[row_ids[static_cast<std::size_t>(
                            row)]];
                    std::uint32_t column_position = descriptor >> 4;
                    const std::uint32_t end = position + (descriptor & 15);
                    const double value = endpoint_values(row);
                    for (; position < end; ++position, ++column_position) {
                        values[m_inner[column_position]] +=
                            m_values[position] * value;
                    }
                }
            };
            if (m_row_patterns32.empty()) {
                apply_patterns(m_row_patterns16);
            } else {
                apply_patterns(m_row_patterns32);
            }
            return;
        }
        const auto apply_rows = [&](const auto& indices, const auto& rows,
                                    auto relative) {
            constexpr bool compact_rows = std::is_same_v<
                typename std::decay_t<decltype(rows)>::value_type,
                std::uint8_t>;
            std::uint32_t position = 0;
            for (Eigen::Index row = 0; row < m_rows; ++row) {
                std::uint32_t begin;
                std::uint32_t end;
                if constexpr (compact_rows) {
                    begin = position;
                    end = begin + rows[static_cast<std::size_t>(row)];
                    position = end;
                } else {
                    begin = rows[static_cast<std::size_t>(row)];
                    end = rows[static_cast<std::size_t>(row) + 1];
                }
                double* row_values = table_values.data();
                if constexpr (decltype(relative)::value) {
                    if (begin != end) {
                        row_values +=
                            m_row_bases[static_cast<std::size_t>(row)];
                    }
                }
                const double value = endpoint_values(row);
                for (std::uint32_t entry = begin; entry < end; ++entry) {
                    row_values[indices[entry]] += m_values[entry] * value;
                }
            }
        };
        if (m_relative_indices) {
            if (m_compact_rows) {
                apply_rows(m_relative_inner, m_row_counts, std::true_type{});
            } else {
                apply_rows(m_relative_inner, m_outer, std::true_type{});
            }
        } else if (m_compact_rows) {
            apply_rows(m_inner, m_row_counts, std::false_type{});
        } else {
            apply_rows(m_inner, m_outer, std::false_type{});
        }
    }

    std::size_t SolarTableInterpolation::storage_bytes() const {
        return m_outer.capacity() * sizeof(std::uint32_t) +
               m_inner.capacity() * sizeof(std::uint32_t) +
               m_row_counts.capacity() * sizeof(std::uint8_t) +
               m_relative_inner.capacity() * sizeof(std::uint16_t) +
               m_row_bases.capacity() * sizeof(std::uint32_t) +
               m_row_patterns16.capacity() * sizeof(std::uint16_t) +
               m_row_patterns32.capacity() * sizeof(std::uint32_t) +
               m_column_patterns.capacity() * sizeof(std::uint32_t) +
               m_values.capacity() * sizeof(double);
    }

    void SolarTransmissionTable::initialize_geometry(
        const std::vector<sasktran2::raytracing::TracedRay>& integration_rays) {
        // find the min/max SZA from the LOS rays and generate the cos_sza_grid
        std::pair<double, double> min_max_cos_sza =
            sasktran2::raytracing::min_max_cos_sza_of_all_rays(
                integration_rays);

        Eigen::VectorXd cos_sza_grid_values;

        if (m_geometry.coordinates().geometry_type() ==
            sasktran2::geometrytype::spherical) {
            // TODO: configure this resolution
            cos_sza_grid_values.setLinSpaced(100, min_max_cos_sza.first,
                                             min_max_cos_sza.second);
        } else {
            // TODO: Can we handle pseudo-spherical here?
            cos_sza_grid_values.resize(1);
            cos_sza_grid_values(0) = min_max_cos_sza.first;
        }

        Eigen::VectorXd alt_values = m_geometry_1d->altitude_grid().grid();

        // create the location interpolator
        m_location_interpolator = std::make_unique<
            sasktran2::grids::AltitudeSZASourceLocationInterpolator>(
            sasktran2::grids::AltitudeGrid(
                std::move(alt_values), sasktran2::grids::gridspacing::constant,
                sasktran2::grids::outofbounds::extend,
                sasktran2::grids::interpolation::linear),
            sasktran2::grids::Grid(std::move(cos_sza_grid_values),
                                   sasktran2::grids::gridspacing::constant,
                                   sasktran2::grids::outofbounds::extend,
                                   sasktran2::grids::interpolation::linear));

        // Construct the matrix that calculates OD on the solar transmission
        // table locations i.e. solar_od_on_grid = matrix @ extinction
        m_geometry_matrix.resize(m_location_interpolator->num_interior_points(),
                                 m_geometry.size());
        m_geometry_matrix.setZero();

        m_ground_hit_flag.resize(
            m_location_interpolator->num_interior_points());

        sasktran2::viewinggeometry::ViewingRay ray_to_sun;

        ray_to_sun.look_away = m_geometry.coordinates().sun_unit();

        raytracing::TracedRay traced_ray;

        for (int i = 0; i < m_location_interpolator->num_interior_points();
             ++i) {
            ray_to_sun.observer.position =
                m_location_interpolator->grid_location(m_geometry.coordinates(),
                                                       i);

            // This method specifically does not allow for refraction
            m_raytracer->trace_ray(ray_to_sun, traced_ray, false);

            if (!traced_ray.ground_is_hit) {
                assign_dense_matrix_column(i, traced_ray, m_geometry_matrix);
                m_ground_hit_flag[i] = false;
            } else {
                m_ground_hit_flag[i] = true;
            }
        }
    }

    void SolarTransmissionTable::generate_interpolation_matrix(
        const std::vector<sasktran2::raytracing::TracedRay>& rays,
        Eigen::SparseMatrix<double, Eigen::RowMajor>& interpolator,
        std::vector<bool>& ground_hit_flag) const {
        // First calculate the number of points we need to create the matrix for
        // We calculate solar transmission at the boundaries of layers, so it is
        // nlayer+1 for each ray
        int numpoints = 0;
        for (const auto& ray : rays) {
            numpoints += (int)ray.layers.size() + 1;
        }

        // od matrix is such that matrix @ extinction = od
        interpolator.resize(numpoints,
                            m_location_interpolator->num_interior_points());

        // Have to handle rays that hit the ground separately since they have no
        // solar transmission
        ground_hit_flag.resize(numpoints, false);

        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;

        std::vector<std::pair<int, double>> interpolator_storage;
        int num_interp;

        int row = 0;
        for (int i = 0; i < rays.size(); ++i) {
            const auto& ray = rays[i];
            for (int j = 0; j < ray.layers.size(); ++j) {
                const auto& layer = ray.layers[j];

                if (j == 0) {
                    // End layer at TOA, need to use layer exit
                    m_location_interpolator->interior_interpolation_weights(
                        m_geometry.coordinates(), layer.exit,
                        interpolator_storage, num_interp);

                    for (int k = 0; k < num_interp; ++k) {
                        tripletList.emplace_back(
                            T(row, interpolator_storage[k].first,
                              interpolator_storage[k].second));
                    }
                    ++row;
                }

                m_location_interpolator->interior_interpolation_weights(
                    m_geometry.coordinates(), layer.entrance,
                    interpolator_storage, num_interp);
                for (int k = 0; k < num_interp; ++k) {
                    tripletList.emplace_back(T(row,
                                               interpolator_storage[k].first,
                                               interpolator_storage[k].second));
                }
                ++row;
            }
        }
        interpolator.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    void SolarTransmissionTable::generate_interpolation(
        const std::vector<sasktran2::raytracing::TracedRay>& rays,
        SolarTableInterpolation& interpolator,
        std::vector<bool>& ground_hit_flag,
        std::vector<Eigen::Vector3d>* solar_propagation_directions) const {
        Eigen::Index numpoints = 0;
        for (const auto& ray : rays) {
            numpoints += static_cast<Eigen::Index>(ray.layers.size()) + 1;
        }
        interpolator.initialize(numpoints,
                                m_location_interpolator->num_interior_points(),
                                numpoints * 4);
        ground_hit_flag.assign(static_cast<std::size_t>(numpoints), false);
        if (solar_propagation_directions != nullptr) {
            solar_propagation_directions->assign(
                static_cast<std::size_t>(numpoints),
                -m_geometry.coordinates().sun_unit());
        }

        std::vector<std::pair<int, double>> weights;
        int num_weights = 0;
        for (const auto& ray : rays) {
            if (ray.layers.empty()) {
                weights.clear();
                interpolator.append_row(weights);
                continue;
            }
            for (int layer_index = 0; layer_index < ray.layers.size();
                 ++layer_index) {
                const auto& layer = ray.layers[layer_index];
                if (layer_index == 0) {
                    m_location_interpolator->interior_interpolation_weights(
                        m_geometry.coordinates(), layer.exit, weights,
                        num_weights);
                    weights.resize(num_weights);
                    interpolator.append_row(weights);
                }
                m_location_interpolator->interior_interpolation_weights(
                    m_geometry.coordinates(), layer.entrance, weights,
                    num_weights);
                weights.resize(num_weights);
                interpolator.append_row(weights);
            }
        }
        interpolator.finalize();
    }

    void
    SolarTransmissionTable::apply(Eigen::Ref<const Eigen::VectorXd> extinction,
                                  Eigen::Ref<Eigen::VectorXd> table_od) const {
        if (extinction.size() != m_geometry_matrix.cols() ||
            table_od.size() != m_geometry_matrix.rows()) {
            throw std::invalid_argument(
                "Invalid 1D solar-table product dimensions");
        }
        table_od.noalias() = m_geometry_matrix * extinction;
    }

    void SolarTransmissionTable::accumulate_transpose(
        Eigen::Ref<const Eigen::VectorXd> table_cotangent,
        Eigen::Ref<Eigen::VectorXd> extinction_cotangent, double scale) const {
        if (table_cotangent.size() != m_geometry_matrix.rows() ||
            extinction_cotangent.size() != m_geometry_matrix.cols()) {
            throw std::invalid_argument(
                "Invalid 1D solar-table transpose dimensions");
        }
        extinction_cotangent.noalias() +=
            scale * m_geometry_matrix.transpose() * table_cotangent;
    }

    std::size_t SolarTransmissionTable::storage_bytes() const {
        return static_cast<std::size_t>(m_geometry_matrix.size()) *
                   sizeof(double) +
               m_ground_hit_flag.capacity() * sizeof(bool);
    }
} // namespace sasktran2::solartransmission
