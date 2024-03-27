#pragma once
#include "sktran_disco/sktran_do.h"

namespace sasktran_disco {

#pragma region "Lazy-Azimuth Framework"
    // A generic object which has some azimuth dependence. Since not all
    // order of the azimuth expansion are needed, this object allows for
    // lazy evaluation of things that depend on the order of the azimuth
    // expansion.
    class AzimuthDependency {
      public:
        // Configures this object for the given azimuth order
        virtual void configureAEOrder(AEOrder m){};

        virtual void postProcessAEOrder(AEOrder m){};
    };

    // We need a way of triggering all azimuth dependent objects from a
    // single call to configureAEOrder (ie. a way of 'cascading' the calls).
    // This object allows dependencies to be registered, when configureAEOrder
    // is called it will call configureAEOrder on all of its registered
    // dependencies.
    class AzimuthDependencyCascade : public AzimuthDependency {
      public:
        // Configures all registered azimuth dependencies of this object for
        // the given azimuth expansion order.
        inline virtual void configureAEOrder(AEOrder m) override {
            for (auto ad : m_dependencies)
                ad->configureAEOrder(m);
        }

        inline virtual void postProcessAEOrder(AEOrder m) override {
            for (auto ad : m_dependencies)
                ad->postProcessAEOrder(m);
        }

        // Register an azimuthal dependency with this object. All subsequent
        // calls to ConfigureAEOrder will also call dependency.ConfigureAEOrder
        inline void registerAzimuthDependency(AzimuthDependency& dependency) {
            m_dependencies.push_back(&dependency);
        }

      private:
        std::list<AzimuthDependency*> m_dependencies;
    };

    // Pure virtual base class which manages the lazy caching of a azimuth
    // dependent function that caches things pre solving the system
    template <class CachedDataType>
    class AzimuthDependentCache : public AzimuthDependency {
      public:
        AzimuthDependentCache() = delete;
        AzimuthDependentCache(uint NSTR)
            : M_NSTR(NSTR), m_cached(M_NSTR, false), m_localdata(M_NSTR),
              m_data(m_localdata), m_reuse_memory_for_all_azimuth(false) {}
        AzimuthDependentCache(uint NSTR, std::vector<CachedDataType>& data)
            : M_NSTR(NSTR), m_cached(M_NSTR, false), m_data(data) {
            if (data.size() == 1) {
                m_reuse_memory_for_all_azimuth = true;
            } else {
                m_reuse_memory_for_all_azimuth = false;
            }
        }

        // Cache the requested azimuth expansion order if it has not already
        // been calculated.
        void configureAEOrder(AEOrder m) override {
            if (m_cached[m] == false) {
                if (m_reuse_memory_for_all_azimuth) {
                    calculateAEOrder(m, m_data[0]);
                } else {
                    calculateAEOrder(m, m_data[m]);
                }
                m_cached[m] = true;
            }
        }

        void reset() {
            for (int i = 0; i < m_cached.size(); ++i) {
                m_cached[i] = false;
            }
        }

        // Cached value accessor.
        inline const CachedDataType& operator[](AEOrder m) const {
            assert(m_cached[m]);
            if (m_reuse_memory_for_all_azimuth) {
                return m_data[0];
            } else {
                return m_data[m];
            }
        }
        // Calculation which is performed when an new order of the cache needs
        // to be calculated.
        virtual void calculateAEOrder(AEOrder m, CachedDataType& val) = 0;

      protected:
        const uint M_NSTR;

      private:
        std::vector<CachedDataType> m_localdata;
        std::vector<CachedDataType>& m_data;
        std::vector<bool> m_cached;

        bool m_reuse_memory_for_all_azimuth;
    };

    // Pure virtual base class which manages the lazy caching of a azimuth
    // dependent
    // function that caches things pre solving the system
    template <class CachedDataType>
    class AzimuthDependentPostCache : public AzimuthDependency {
      public:
        AzimuthDependentPostCache() = delete;
        AzimuthDependentPostCache(uint NSTR)
            : M_NSTR(NSTR), m_cached(M_NSTR, false), m_data(M_NSTR) {}
        // Cache the requested azimuth expansion order if it has not already
        // been calculated.
        void postProcessAEOrder(AEOrder m) override {
            if (m_cached[m] == false) {
                calculateAEOrder(m, m_data[m]);
                m_cached[m] = true;
            }
        }
        // Cached value accessor.
        inline const CachedDataType& operator[](AEOrder m) const {
            assert(m_cached[m]);
            return m_data[m];
        }
        // Calculation which is performed when an new order of the cache needs
        // to be calculated.
        virtual void calculateAEOrder(AEOrder m, CachedDataType& val) = 0;

      protected:
        const uint M_NSTR;

      private:
        std::vector<CachedDataType> m_data;
        std::vector<bool> m_cached;
    };

#pragma endregion

#pragma region "Lazy-Azimuth Cache Objects"
    class AlbedoExpansion : public AzimuthDependentCache<Albedo> {
      public:
        AlbedoExpansion(const std::vector<LineOfSight>& los,
                        const std::vector<double>& streams, double csz,
                        std::unique_ptr<BRDF_Base> brdf, uint nterms)
            : AzimuthDependentCache<Albedo>(static_cast<uint>(streams.size())),
              M_LOS(los), M_MU(streams), M_CSZ(csz), m_brdf(std::move(brdf)),
              m_nterms(nterms) {
            // empty
        }

        void injectTestingBRDF(std::unique_ptr<BRDF_Base> brdf) {
            m_brdf = std::move(brdf);
        }

        virtual void calculateAEOrder(AEOrder m,
                                      Albedo& sum_matrix) override final;

      private:
        const std::vector<LineOfSight>& M_LOS;
        const std::vector<double>& M_MU;
        const double M_CSZ;
        std::unique_ptr<BRDF_Base> m_brdf;
        uint m_nterms;
    };

} // namespace sasktran_disco
#pragma endregion
