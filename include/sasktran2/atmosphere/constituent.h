#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/grids.h>
#include <sasktran2/geometry.h>
#include <sasktran2/atmosphere/grid_storage.h>

namespace sasktran2::atmosphere {
    class Constituent {
      private:
        std::unique_ptr<AtmosphereGridStorage> m_storage;
        Eigen::VectorXd m_location_scale_factor;

        virtual void
        internal_populate_stokes1(AtmosphereGridStorageFull<1>& storage,
                                  const Geometry& geometry) = 0;
        virtual void
        internal_populate_stokes3(AtmosphereGridStorageFull<3>& storage,
                                  const Geometry& geometry) = 0;

      public:
        void populate_storage(const Geometry& geometry, int nstokes) const;

        template <int NSTOKES>
        void add_to_storage(AtmosphereGridStorageFull<NSTOKES>& storage);

        const std::unique_ptr<AtmosphereGridStorage>& storage() const {
            return m_storage;
        }
    };

    template <int NSTOKES> class ManualConstituent : public Constituent {
      private:
        virtual void
        internal_populate_stokes1(AtmosphereGridStorageFull<1>& storage,
                                  const Geometry& geometry){};
        virtual void
        internal_populate_stokes3(AtmosphereGridStorageFull<3>& storage,
                                  const Geometry& geometry){};

      public:
    };

} // namespace sasktran2::atmosphere
