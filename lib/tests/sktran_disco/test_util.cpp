#include "sasktran2/atmosphere/atmosphere.h"
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
     * lowlevel interface instead.
     *
     * @tparam NSTOKES
     * @param test_case
     * @param test_spec
     * @param radiance
     * @param wf
     */
    template <int NSTOKES>
    void run_lowlevel_from_old_testspec(
        sasktran_disco::testing::TestCase<NSTOKES>& test_case,
        sasktran_disco::SKTRAN_DO_TestSpec<NSTOKES, -1>* test_spec,
        std::vector<double>& radiance, std::vector<double>& abs_diff,
        std::vector<std::vector<double>>* wf) {
        // Set up the low level interface
        int nstr = test_case.nstr;
        int nlyr = test_case.nlyr;
        int nlos = (int)test_case.linesofsight.size();
        int nderiv;

        if (test_spec) {
            int nderiv = (int)test_spec->perturbations()->size();
        } else {
            nderiv = 0;
        }

        sasktran_disco_lowlevel::CPPApi cppapi(nstr, nlyr, 1, NSTOKES, nlos,
                                               nderiv);

        // Copy the atmosphere
        for (int l = 0; l < nlyr; ++l) {
            cppapi.od(l, 0) = test_case.layers[l].optical_depth;
            cppapi.ssa(l, 0) = test_case.layers[l].ssa;

            for (int k = 0; k < nstr; ++k) {
                cppapi.a1(k, l, 0) = test_case.layers[l].lephasef[k].a1;

                if constexpr (NSTOKES == 3) {
                    cppapi.a2(k, l, 0) = test_case.layers[l].lephasef[k].a2;
                    cppapi.a3(k, l, 0) = test_case.layers[l].lephasef[k].a3;
                    cppapi.b1(k, l, 0) = test_case.layers[l].lephasef[k].b1;
                }
            }
        }

        // Now the albedo
        cppapi.albedo(0) = test_case.lambertian;

        // Copy the solar zenith
        cppapi.cos_sza() = test_case.solar.csz;

        // Set the layer boundaries, not used since we aren't using pseudo
        // spherical anyways
        for (int l = 0; l < nlyr; ++l) {
            cppapi.layerboundaryaltitude(l) =
                100000.0 - double((100000.0) * l) / double(nlyr);
        }
        cppapi.earth_radius() = 6372000;

        // And the LOS angles
        for (int i = 0; i < nlos; ++i) {
            cppapi.cos_vza(i) = test_case.linesofsight[i].coszen;
            cppapi.saa(i) = test_case.solar.saz - test_case.linesofsight[i].az;
        }

        sasktran_disco_lowlevel::Config config;
        sasktran_disco_lowlevel::Output output;
        sasktran_disco_lowlevel::Atmosphere atmosphere;
        sasktran_disco_lowlevel::ViewingGeometry viewinggeo;
        sasktran_disco_lowlevel::WeightingFunctions weightingfunctions;

        cppapi.initialize_c_api(&atmosphere, &config, &viewinggeo,
                                &weightingfunctions, &output);

        // Set additional config options
        if (test_spec) {
            config.numazimuthexpansion =
                test_spec->getForcedNumberAzimuthTerms();
        }
        config.useexactsinglescatter = false;
        config.usepseudospherical = false;

        sasktran_disco_lowlevel::calculate(
            &atmosphere, &config, &weightingfunctions, &viewinggeo, &output);

        // Copy the results to the old output format
        radiance.resize(nlos * NSTOKES);
        abs_diff.resize(nlos * NSTOKES);
        for (int i = 0; i < nlos * NSTOKES; ++i) {
            radiance[i] =
                output.radiance[i] * test_case.solar.intensities.direct;
            abs_diff[i] = abs(radiance[i] - (*test_case.correct_radiances)[i]);
        }
    }

    template void run_lowlevel_from_old_testspec(
        sasktran_disco::testing::TestCase<1>& test_case,
        sasktran_disco::SKTRAN_DO_TestSpec<1, -1>* test_spec,
        std::vector<double>& radiance, std::vector<double>& abs_diff,
        std::vector<std::vector<double>>* wf);

    template void run_lowlevel_from_old_testspec(
        sasktran_disco::testing::TestCase<3>& test_case,
        sasktran_disco::SKTRAN_DO_TestSpec<3, -1>* test_spec,
        std::vector<double>& radiance, std::vector<double>& abs_diff,
        std::vector<std::vector<double>>* wf);

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

        if (test_spec) {
            config.set_num_do_forced_azimuth(
                test_spec->getForcedNumberAzimuthTerms());
        }

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
        sasktran2::atmosphere::Surface surface;

        surface.albedo().resize(nwavel);
        surface.albedo().setConstant(test_case.lambertian);

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
