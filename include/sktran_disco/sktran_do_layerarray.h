#pragma once
#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_surface.h"
#include "sktran_disco/sktran_do_types.h"
#include <sasktran2/atmosphere/atmosphere.h>

namespace sasktran_disco {
    template <Propagating dir, int NSTOKES, int CNSTR = -1>
    class OpticalLayerArrayIterator;

    template <int NSTOKES, int CNSTR = -1>
    using OpticalLayerArrayROP =
        ReadOnlyProperties<BasicProperties<NSTOKES>, SolarProperties<NSTOKES>,
                           UserSpecProperties>;

    /**
     * @brief A vector of OpticalLayer objects which represent the layers of the
     * atmosphere.
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR = -1>
    class OpticalLayerArray : public OpticalLayerArrayROP<NSTOKES> {
      public:
        /**
         * @brief Construct a new Optical Layer Array object using the sasktran2
         * interface This is the primary way of constructing the
         * OpticalLayerArray.
         *
         * @param config
         * @param wavelidx
         * @param los
         * @param brdf
         * @param geometry_layers
         * @param atmosphere
         * @param sk_config
         */
        OpticalLayerArray(
            const PersistentConfiguration<NSTOKES, CNSTR>& config, int wavelidx,
            const std::vector<LineOfSight>& los,
            const GeometryLayerArray<NSTOKES, CNSTR>& geometry_layers,
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
            const sasktran2::Config& sk_config);

        /**
         * @brief Internal method which is called to compute the reflected
         * intensities for the line of sight rays.
         *
         * @param m
         * @param los
         */
        void
        computeReflectedIntensities(AEOrder m,
                                    const sasktran_disco::LineOfSight& los);

        /**
         * @brief Layer at lidx
         *
         * @param lidx
         * @return OpticalLayer<NSTOKES, CNSTR>&
         */
        inline OpticalLayer<NSTOKES, CNSTR>& operator[](uint lidx) {
            return *m_layers[lidx];
        }

        /**
         * @brief Layer at lidx
         *
         * @param lidx
         * @return const OpticalLayer<NSTOKES, CNSTR>&
         */
        inline const OpticalLayer<NSTOKES, CNSTR>& operator[](uint lidx) const {
            return *m_layers[lidx];
        }

        /**
         * @brief Bottom layer of the atmosphere
         *
         * @return OpticalLayer<NSTOKES, CNSTR>&
         */
        inline OpticalLayer<NSTOKES, CNSTR>& bottom() {
            return *m_layers.back();
        }

        /**
         * @brief Bottom layer of the atmosphere
         *
         * @return const OpticalLayer<NSTOKES, CNSTR>&
         */
        inline const OpticalLayer<NSTOKES, CNSTR>& bottom() const {
            return *m_layers.back();
        }

        /**
         * @brief Layer at lidx
         *
         * @param lidx
         * @return OpticalLayer<NSTOKES, CNSTR>&
         */
        inline OpticalLayer<NSTOKES, CNSTR>& layer(uint lidx) {
            return *m_layers[lidx];
        }

        /**
         * @brief Layer at lidx
         *
         * @param lidx
         * @return const OpticalLayer<NSTOKES, CNSTR>&
         */
        inline const OpticalLayer<NSTOKES, CNSTR>& layer(uint lidx) const {
            return *m_layers[lidx];
        }

        /**
         * @brief Layer at the top of the atmosphere
         *
         * @return OpticalLayer<NSTOKES, CNSTR>&
         */
        inline OpticalLayer<NSTOKES, CNSTR>& top() { return *m_layers[0]; }

        /**
         * @brief Layer at the top of the atmosphere
         *
         * @return const OpticalLayer<NSTOKES, CNSTR>&
         */
        inline const OpticalLayer<NSTOKES, CNSTR>& top() const {
            return *m_layers[0];
        }

        /**
         * @brief The solar irradiance at the top of the atmosphere
         *
         * @return double
         */
        inline double directIntensityTOA() const { return m_direct_toa; }

        /**
         * @brief The number of layers in the atmosphere
         *
         * @return uint
         */
        inline uint numLayers() const { return this->M_NLYR; }

        /**
         * @brief The user input derivatives
         *
         * @return const InputDerivatives<NSTOKES>&
         */
        inline const InputDerivatives<NSTOKES>& inputDerivatives() const {
            return m_input_derivatives;
        }

        /**
         * @brief The user input derivatives
         *
         * @return InputDerivatives<NSTOKES>&
         */
        inline InputDerivatives<NSTOKES>& inputDerivatives() {
            return m_input_derivatives;
        }

        /**
         * @brief The surface
         *
         * @return Surface<NSTOKES, CNSTR>&
         */
        inline Surface<NSTOKES, CNSTR>& surface() { return m_surface; }

        /**
         * @brief The total reflected intensity for a given line of sight
         *     and azimuth order
         *
         * @param m
         * @param los
         * @return const Radiance<NSTOKES>&
         */
        inline const Radiance<NSTOKES>&
        reflectedIntensity(AEOrder m, const LineOfSight& los) {
            if (!m_reflection_computed[m][los.unsorted_index]) {
                computeReflectedIntensities(m, los);
            }
            return m_ground_reflection[m][los.unsorted_index];
        }

        /**
         * @brief Returns an iterator which will iterate in the specified
         * direction to the given terminal optical depth.
         *
         * @tparam dir
         * @param terminal
         * @return OpticalLayerArrayIterator<dir, NSTOKES, CNSTR>
         */
        template <Propagating dir>
        inline OpticalLayerArrayIterator<dir, NSTOKES, CNSTR>
        iteratorTo(double terminal) {
            return OpticalLayerArrayIterator<dir, NSTOKES, CNSTR>(*this,
                                                                  terminal);
        }

        /**
         * @brief Returns an iterator which will go through all layers in the
         * specified direction
         *
         * @tparam dir
         * @return OpticalLayerArrayIterator<dir, NSTOKES, CNSTR>
         */
        template <Propagating dir>
        inline OpticalLayerArrayIterator<dir, NSTOKES, CNSTR> iteratorAcross() {
            if (dir == Propagating::UP) {
                return iteratorTo<dir>(0);
            } else {
                return iteratorTo<dir>(bottom().opticalDepth(Location::FLOOR));
            }
        }

        /**
         * @brief Returns the first layer which has a total optical depth
         * greater than the given value
         *
         * @param optical_depth
         * @return OpticalLayer<NSTOKES, CNSTR>*
         */
        inline OpticalLayer<NSTOKES, CNSTR>* layerAt(double optical_depth) {
            for (LayerIndex p = 0; p < this->M_NLYR; ++p) {
                if (layer(p).opticalDepth(Location::FLOOR) >= optical_depth)
                    return &layer(p);
            }
            return nullptr;
        }

        /**
         * @brief Returns the first layer which has a total optical depth
         * greater than the given value
         *
         * @param optical_depth
         * @return const OpticalLayer<NSTOKES, CNSTR>*
         */
        inline const OpticalLayer<NSTOKES, CNSTR>*
        layerAt(double optical_depth) const {
            for (LayerIndex p = 0; p < this->M_NLYR; ++p) {
                if (layer(p).opticalDepth(Location::FLOOR) >= optical_depth)
                    return &layer(p);
            }
            return nullptr;
        }

        /**
         * @brief Returns the layer which contains the altitude given
         *
         * @param altitude
         * @return const OpticalLayer<NSTOKES, CNSTR>*
         */
        inline const OpticalLayer<NSTOKES, CNSTR>*
        layerAtAltitude(double altitude) const {
            // Binary search since this can be slow and called many times

            int lower = 0;
            int upper = this->M_NLYR - 1;

            int retindex = 0;
            while (true) {
                if (lower == upper) {
                    retindex = upper;
                    break;
                }
                if (upper - lower == 1) {
                    retindex =
                        layer(lower).altitude(Location::FLOOR) <= altitude
                            ? lower
                            : upper;
                    break;
                }
                int checkIndex = (lower + upper) / 2;

                if (layer(checkIndex).altitude(Location::FLOOR) > altitude) {
                    lower = checkIndex;
                } else {
                    upper = checkIndex;
                }
            }
            return &layer(retindex);
        }

        /**
         * @brief Calculates the optical depth at the given altitude
         *
         * @param altitude
         * @return double
         */
        inline double opticalDepthAt(double altitude) const {
            const auto& layer = layerAtAltitude(altitude);

            if (altitude > layer->altitude(Location::CEILING)) {
                return 0.0;
            } else {
                double f = 1 - (layer->altitude(Location::CEILING) - altitude) /
                                   (layer->altitude(Location::CEILING) -
                                    layer->altitude(Location::FLOOR));

                return layer->opticalDepth(Location::FLOOR) -
                       layer->opticalDepth(Location::INSIDE) * f;
            }
        }

        /**
         * @brief Returns the wavelength index that the layer array is loaded
         * with
         *
         * @return int
         */
        int wavelength_index() const { return m_wavel_index; }

      protected:
        /**
         * @brief Calculates the transmission factors for the layers
         *
         */
        void configureTransmission();

      protected:
        VectorDim1<std::unique_ptr<OpticalLayer<NSTOKES, CNSTR>>>
            m_layers;        /** Layers */
        double m_direct_toa; /** Direct solar TOA irradiance */
        size_t m_num_los;    /** Number of lines of sight */

        InputDerivatives<NSTOKES>&
            m_input_derivatives; /** User input derivatives */
        size_t m_wavel_index; /** Index to the loaded wavelength when using the
                                 sasktran2 interface */

        VectorDim2<Radiance<NSTOKES>>
            m_ground_reflection; /** Vector of ground reflection objects
                                    [azi][los]*/
        VectorDim2<bool>
            m_reflection_computed; /** True if the ground reflection has been
                                      computed [azi][los]*/

        const Eigen::MatrixXd&
            m_chapman_factors; /** Chapman factors for the layers*/

        const Eigen::MatrixXd&
            m_optical_interpolator; /** Interpolation matrix from the atmosphere
                                       levels to the layers */

        std::vector<Dual<double>>& m_transmission; /** Transmission factors for
                                                     the layers */

        const PersistentConfiguration<NSTOKES, CNSTR>&
            m_config; /** Internal config object */

        Surface<NSTOKES, CNSTR> m_surface;

        bool m_include_direct_bounce; /** True if when computing the reflected
                                         intensities the direct solar bounce
                                         should be included */
    };

    // STL random access iterator with additional layer of functionality for
    // iterating over a range of a OpticalLayerArray.
    template <Propagating dir, int NSTOKES, int CNSTR>
    class OpticalLayerArrayIterator {
      public:
        // Constructors
        OpticalLayerArrayIterator() noexcept
            : m_layers(nullptr), m_current_layer_index(-1),
              m_terminal_optical_depth(std::nan("1")) {
            // empty
        }
        OpticalLayerArrayIterator(
            const OpticalLayerArrayIterator& other) noexcept = default;
        OpticalLayerArrayIterator(OpticalLayerArrayIterator&& other) noexcept =
            default;
        OpticalLayerArrayIterator&
        operator=(OpticalLayerArrayIterator&& other) noexcept = default;
        OpticalLayerArrayIterator&
        operator=(const OpticalLayerArrayIterator& other) noexcept = default;

        // STL type info
        typedef int difference_type;
        typedef OpticalLayer<NSTOKES, CNSTR> value_type;
        typedef OpticalLayer<NSTOKES, CNSTR>* pointer;
        typedef OpticalLayer<NSTOKES, CNSTR>& reference;
        typedef std::random_access_iterator_tag iterator_category;

        // Comparison operators
        bool operator==(const OpticalLayerArrayIterator& other) const {
            return m_layers == other.m_layers &&
                   m_current_layer_index == other.m_current_layer_index &&
                   m_terminal_optical_depth == other.m_terminal_optical_depth;
        }
        inline bool operator!=(const OpticalLayerArrayIterator& other) const {
            return *this == other;
        }
        bool operator<(const OpticalLayerArrayIterator& other) const {
            if (isDownwelling()) {
                return m_current_layer_index < other.m_current_layer_index;
            } else {
                return m_current_layer_index > other.m_current_layer_index;
            }
        }
        bool operator<=(const OpticalLayerArrayIterator& other) const {
            return *this < other || *this == other;
        }
        bool operator>(const OpticalLayerArrayIterator& other) const {
            return !(*this <= other);
        }
        bool operator>=(const OpticalLayerArrayIterator& other) const {
            return !(*this < other);
        }

        // Dereference operators
        OpticalLayer<NSTOKES, CNSTR>&
        operator*() { // needs to be removed for const
            return layer();
        }
        const OpticalLayer<NSTOKES, CNSTR>& operator*() const {
            return layer();
        }
        OpticalLayer<NSTOKES, CNSTR>&
        operator->() { // needs to be removed for const
            return layer();
        }
        const OpticalLayer<NSTOKES, CNSTR>& operator->() const {
            return layer();
        }
        inline OpticalLayer<NSTOKES, CNSTR>& operator[](int offset) {
            return m_layers->layer(m_current_layer_index + offset);
        }
        inline const OpticalLayer<NSTOKES, CNSTR>&
        operator[](int offset) const {
            return m_layers->layer(m_current_layer_index + offset);
        }

        // Increment/decrement operators
        inline OpticalLayerArrayIterator& operator++() {
            if (isDownwelling())
                ++m_current_layer_index;
            else
                --m_current_layer_index;
            return *this;
        }
        inline OpticalLayerArrayIterator& operator++(int) {
            if (isDownwelling())
                ++m_current_layer_index;
            else
                --m_current_layer_index;
            return *this;
        }
        inline OpticalLayerArrayIterator& operator--() {
            if (isDownwelling())
                --m_current_layer_index;
            else
                ++m_current_layer_index;
            return *this;
        }
        inline OpticalLayerArrayIterator& operator--(int) {
            if (isDownwelling())
                --m_current_layer_index;
            else
                ++m_current_layer_index;
            return *this;
        }
        inline OpticalLayerArrayIterator& operator+=(int offset) {
            if (isDownwelling())
                m_current_layer_index += offset;
            else
                m_current_layer_index -= offset;
            return *this;
        }
        inline OpticalLayerArrayIterator& operator-=(int offset) {
            return operator+=(-offset);
        }

        // Arithmetic operators
        inline OpticalLayerArrayIterator operator+(int offset) const {
            return OpticalLayerArrayIterator(*this) += offset;
        }
        inline OpticalLayerArrayIterator operator-(int offset) const {
            return OpticalLayerArrayIterator(*this) -= offset;
        }
        inline int operator+(const OpticalLayerArrayIterator& other) const {
            return m_current_layer_index + other.m_current_layer_index;
        }
        inline int operator-(const OpticalLayerArrayIterator& other) const {
            return m_current_layer_index - other.m_current_layer_index;
        }

      public: // Interface
        // Constructs a iterator over all layers.
        OpticalLayerArrayIterator(OpticalLayerArray<NSTOKES, CNSTR>& layers) {
            m_layers = &layers;
            if (isDownwelling()) {
                m_current_layer_index = 0;
            } else {
                m_current_layer_index = layers.numLayers() - 1;
            }
            if (isDownwelling()) {
                m_terminal_optical_depth =
                    layers.bottom().opticalDepth(Location::FLOOR);
            } else {
                m_terminal_optical_depth =
                    layers.top().opticalDepth(Location::CEILING);
            }
        }
        // Constructs an iterator which terminates at terminal_depth.
        OpticalLayerArrayIterator(OpticalLayerArray<NSTOKES, CNSTR>& layers,
                                  double terminal_depth)
            : OpticalLayerArrayIterator(layers) {
            m_terminal_optical_depth = terminal_depth;
        }

        // Returns non-mutable OpticalLayer reference
        inline const OpticalLayer<NSTOKES, CNSTR>& layer() const {
            return m_layers->layer(m_current_layer_index);
        }
        // Returns mutable OpticalLayer reference
        inline OpticalLayer<NSTOKES, CNSTR>& layer() {
            return m_layers->layer(m_current_layer_index);
        }
        // Returns true if this iterator is still accessing a valid range of the
        // OpticalLayerArray false otherwise. It is important to not perform any
        // operations an an iterator when isValid() is false. This indicates
        // that the iterator has expired.
        inline bool isValid() const {
            bool ok = true;
            ok &= m_current_layer_index <
                      static_cast<int>(m_layers->numLayers()) &&
                  m_current_layer_index >= 0;
            if (!ok)
                return ok;
            if (isUpwelling()) {
                ok &= layer().opticalDepth(Location::FLOOR) >
                      m_terminal_optical_depth;
            } else {
                ok &= layer().opticalDepth(Location::CEILING) <
                      m_terminal_optical_depth;
            }
            return ok;
        }

        // Return the entry point to a layer (upwelling is FLOOR, down welling
        // is CEILING)
        inline Location entry() const {
            return isUpwelling() ? Location::FLOOR : Location::CEILING;
        }
        // Returns the exit point to a layer (see entry())
        inline Location exit() const {
            if (m_terminal_optical_depth >
                    layer().opticalDepth(Location::FLOOR) &&
                m_terminal_optical_depth <
                    layer().opticalDepth(Location::CEILING)) {
                return Location::INSIDE;
            } else {
                return isUpwelling() ? Location::CEILING : Location::FLOOR;
            }
        }
        // Returns the optical depth at entry to the current layer
        inline double entryOpticalDepth() const {
            return isUpwelling() ? layer().opticalDepth(Location::FLOOR)
                                 : layer().opticalDepth(Location::CEILING);
        }
        // Returns the optical depth at the exit of the current layer
        inline double exitOpticalDepth() const {
            if (exit() == Location::INSIDE)
                return m_terminal_optical_depth;
            else {
                return isUpwelling() ? layer().opticalDepth(Location::CEILING)
                                     : layer().opticalDepth(Location::FLOOR);
            }
        }
        // Returns true if the iterator is upwelling false otherwise.
        constexpr inline bool isUpwelling() const {
            return dir == Propagating::UP;
        }
        // Returns true if the iterator is down welling false otherwise.
        constexpr inline bool isDownwelling() const {
            return dir == Propagating::DOWN;
        }

        // Returns mutable pointer to the layer
        inline OpticalLayer<NSTOKES, CNSTR>* ptr() { return &layer(); }

        // Returns non-mutable pointer to the layer
        inline const OpticalLayer<NSTOKES, CNSTR>* ptr() const {
            return &layer();
        }

      private:
        double m_terminal_optical_depth;
        int m_current_layer_index;
        OpticalLayerArray<NSTOKES, CNSTR>* m_layers;
    };
} // namespace sasktran_disco
