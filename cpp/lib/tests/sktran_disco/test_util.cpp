#include "sasktran2/atmosphere/atmosphere.h"
#include "sasktran2/atmosphere/surface.h"
#include "sasktran2/config.h"
#include "sasktran2/geometry.h"
#include "sasktran2/grids.h"
#include "sasktran2/viewinggeometry.h"
#include <sasktran2/test_helper.h>

#include <sasktran2/testing/do_test_util.h>

#include <sasktran2.h>

namespace sasktran_disco::testing {
    /** This is an interface function to run the old SASKTRAN DISCO tests that
     * used the sasktran_disco::testing::TestCase class along with the standard
     * DO engine.  Here we transform the TestCase input to input used in the
     * sasktran2 interface instead
     *
     * @tparam NSTOKES
     * @param test_case
     * @param test_spec
     * @param radiance
     * @param wf
     */
    template <int NSTOKES>
    void run_sasktran2_from_old_testspec(
        sasktran_disco::testing::TestCase<NSTOKES>& test_case,
        sasktran_disco::SKTRAN_DO_TestSpec<NSTOKES, -1>* test_spec,
        std::vector<double>& radiance, std::vector<double>& abs_diff,
        std::vector<std::vector<double>>* wf) {
        // Set up the low level interface
        int nlyr = test_case.nlyr;
        int nlos = (int)test_case.linesofsight.size();

        auto config = sasktran2::Config();
        config.set_num_do_streams(test_case.nstr);
        config.set_num_stokes(NSTOKES);
        config.set_multiple_scatter_source(
            sasktran2::Config::MultipleScatterSource::discrete_ordinates);
        config.set_single_scatter_source(
            sasktran2::Config::SingleScatterSource::discrete_ordinates);

        Eigen::VectorXd grid_values(nlyr + 1);
        for (int i = 0; i < nlyr + 1; ++i) {
            grid_values(i) = i;
        }

        auto viewing_geo =
            sasktran2::viewinggeometry::ViewingGeometryContainer();

        for (int i = 0; i < nlos; ++i) {
            viewing_geo.add_ray(sasktran2::viewinggeometry::GroundViewingSolar(
                test_case.solar.csz,
                -test_case.solar.saz + test_case.linesofsight[i].az,
                test_case.linesofsight[i].coszen,
                grid_values(Eigen::last) + 1));
        }

        // Construct the Atmosphere
        int nwavel = 1;
        sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES> storage(
            nwavel, grid_values.size(), test_case.nstr);
        sasktran2::atmosphere::Surface<NSTOKES> surface(nwavel);

        surface.brdf_args().setConstant(test_case.lambertian);

        auto geometry = sasktran2::Geometry1D(
            test_case.solar.csz, 0, 6372000, std::move(grid_values),
            sasktran2::grids::interpolation::lower,
            sasktran2::geometrytype::planeparallel);

        // Copy the atmosphere
        for (int l = 0; l < nlyr; ++l) {
            int atmidx = nlyr - l - 1;

            storage.total_extinction(atmidx, 0) =
                test_case.layers[l].optical_depth;
            storage.ssa(atmidx, 0) = test_case.layers[l].ssa;

            for (int k = 0; k < test_case.nstr; ++k) {

                if constexpr (NSTOKES == 3) {
                    storage.leg_coeff(k * 4, atmidx, 0) =
                        test_case.layers[l].lephasef[k].a1;
                    storage.leg_coeff(k * 4 + 1, atmidx, 0) =
                        test_case.layers[l].lephasef[k].a2;
                    storage.leg_coeff(k * 4 + 2, atmidx, 0) =
                        test_case.layers[l].lephasef[k].a3;
                    storage.leg_coeff(k * 4 + 3, atmidx, 0) =
                        -test_case.layers[l].lephasef[k].b1;
                } else {
                    storage.leg_coeff(k, atmidx, 0) =
                        test_case.layers[l].lephasef[k].a1;
                }
            }
        }
        sasktran2::atmosphere::Atmosphere<NSTOKES> atmo(
            std::move(storage), std::move(surface), true);

        // Make the engine
        Sasktran2<NSTOKES> engine(config, &geometry, viewing_geo);

        sasktran2::OutputIdealDense<NSTOKES> output;
        engine.calculate_radiance(atmo, output);

        // Copy the results to the old output format
        radiance.resize(nlos * NSTOKES);
        abs_diff.resize(nlos * NSTOKES);
        for (int i = 0; i < nlos * NSTOKES; ++i) {
            radiance[i] =
                output.radiance().value(i) * test_case.solar.intensities.direct;
            abs_diff[i] = abs(radiance[i] - (*test_case.correct_radiances)[i]);
        }
    }

    template void run_sasktran2_from_old_testspec(
        sasktran_disco::testing::TestCase<1>& test_case,
        sasktran_disco::SKTRAN_DO_TestSpec<1, -1>* test_spec,
        std::vector<double>& radiance, std::vector<double>& abs_diff,
        std::vector<std::vector<double>>* wf);

    template void run_sasktran2_from_old_testspec(
        sasktran_disco::testing::TestCase<3>& test_case,
        sasktran_disco::SKTRAN_DO_TestSpec<3, -1>* test_spec,
        std::vector<double>& radiance, std::vector<double>& abs_diff,
        std::vector<std::vector<double>>* wf);

} // namespace sasktran_disco::testing
