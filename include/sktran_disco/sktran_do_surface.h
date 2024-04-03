#pragma once
#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_types.h"
#include <sasktran2/atmosphere/surface.h>

namespace sasktran_disco {
    template <int NSTOKES, int CNSTR = -1> struct SurfaceStorage {
        static constexpr int CSIZE = CNSTR == -1 ? Eigen::Dynamic : CNSTR;

        struct SurfaceExpansion {
            Eigen::Matrix<double, -1, CSIZE> los_stream;
            Eigen::Matrix<double, CSIZE, CSIZE> stream_stream;
            Eigen::Matrix<double, -1, 1> los_solar;
            Eigen::Matrix<double, CSIZE, 1> stream_solar;
        };

        SurfaceExpansion brdf;
        std::vector<SurfaceExpansion> d_brdf;

        Eigen::VectorXd m_quadrature_phi;
        Eigen::VectorXd m_quadrature_weight;

        // Geometry factors
        double csz;                     /** Cosine solar zenith */
        const std::vector<double>* mu;  /** Stream angles */
        const std::vector<double>* mu0; /** LOS angles */

        void resize(int nstr, int nlos, int nderiv) {
            brdf.los_stream.resize(nlos, nstr / 2);
            brdf.stream_stream.resize(nstr / 2, nstr / 2);
            brdf.los_solar.resize(nlos, 1);
            brdf.stream_solar.resize(nstr / 2, 1);

            d_brdf.resize(nderiv);
            for (auto& d_b : d_brdf) {
                d_b.los_stream.resize(nlos, nstr / 2);
                d_b.stream_stream.resize(nstr / 2, nstr / 2);
                d_b.los_solar.resize(nlos, 1);
                d_b.stream_solar.resize(nstr / 2, 1);
            }
        }
    };

    template <int NSTOKES, int CNSTR = -1> class Surface {
      private:
        SurfaceStorage<NSTOKES, CNSTR>& m_storage;

        const sasktran2::atmosphere::Surface<NSTOKES>& m_sk2_surface;
        int m_wavel_idx;

      public:
        Surface(SurfaceStorage<NSTOKES, CNSTR>& storage,
                const sasktran2::atmosphere::Surface<NSTOKES>& sk2_surface,
                int wavel_idx)
            : m_storage(storage), m_sk2_surface(sk2_surface),
              m_wavel_idx(wavel_idx) {}

        const sasktran2::atmosphere::Surface<NSTOKES>& sk2_surface() const {
            return m_sk2_surface;
        }

        void calculate(AEOrder m) {
            double albedo =
                m_sk2_surface.brdf(m_wavel_idx, 0, 0, 0)(0, 0) * EIGEN_PI;

            m_storage.brdf.los_stream.setConstant(albedo);
            m_storage.brdf.stream_solar.setConstant(albedo);
            m_storage.brdf.stream_stream.setConstant(albedo);
            m_storage.brdf.los_solar.setConstant(albedo);
        }

        const SurfaceStorage<NSTOKES, CNSTR>& storage() const {
            return m_storage;
        }
        SurfaceStorage<NSTOKES, CNSTR>& storage() { return m_storage; }
    };

} // namespace sasktran_disco
