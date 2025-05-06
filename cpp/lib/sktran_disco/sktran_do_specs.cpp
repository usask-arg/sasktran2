#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_specs.h"
#include "sktran_disco/sktran_do_quadrature.h"
#include <sasktran2/math/wigner.h>

sasktran_disco::SKTRAN_DO_UserSpec::SKTRAN_DO_UserSpec() {
    configureDefaultDetails();
}

void sasktran_disco::SKTRAN_DO_UserSpec::configureDefaultDetails() {
    m_nstr = 0;
    m_nlyr = 0;
    setSSAEqual1Dither();
    setTOAIntensities();
    setSSOnly();
}

#pragma region "Setters"

void sasktran_disco::SKTRAN_DO_UserSpec::configure(
    unsigned int num_streams, unsigned int num_atmo_layers,
    double solar_direct_intensity) {
    setNumberOfStreams(num_streams);
    cacheLPOfStreamAngles();
    setNumberOfLayers(num_atmo_layers);
    setTOAIntensities(solar_direct_intensity);
}

sasktran_disco::SKTRAN_DO_UserSpec*
sasktran_disco::SKTRAN_DO_UserSpec::setNumberOfStreams(
    unsigned int num_streams) {
    using namespace sasktran_disco;
    // check that number of streams is GE 4 and even
    if (num_streams < 2 || num_streams % 2 == 1 || num_streams > 40) {
        throw InvalidConfiguration(
            "Number of streams must be an even number in the range [2, 40]!");
    }
    m_nstr = (uint)num_streams;
    return this;
}

sasktran_disco::SKTRAN_DO_UserSpec*
sasktran_disco::SKTRAN_DO_UserSpec::setSSOnly(bool ss_only) {
    m_ss_only = ss_only;

    return this;
}

void sasktran_disco::SKTRAN_DO_UserSpec::cacheLPOfStreamAngles() {
    using namespace sasktran_disco;
    // get stream angles ang weights from tables
    getStreamsAndWeights(m_nstr, m_abscissae, m_weights);
    // cache legendre polynomials evaluated at stream angles
    m_lp_abscissae4.resize(
        m_nstr,
        VectorDim2<sasktran_disco::LegendrePhaseContainer<4>>(
            m_nstr,
            VectorDim1<sasktran_disco::LegendrePhaseContainer<4>>(m_nstr)));
    m_lp_abscissae1.resize(
        m_nstr,
        VectorDim2<sasktran_disco::LegendrePhaseContainer<1>>(
            m_nstr,
            VectorDim1<sasktran_disco::LegendrePhaseContainer<1>>(m_nstr)));

    for (AEOrder m = 0; m < m_nstr; ++m) {
        auto calculatorP = sasktran2::math::WignerDCalculator(m, 0);
        auto calculatorneg = sasktran2::math::WignerDCalculator(m, -2);
        auto calculatorpos = sasktran2::math::WignerDCalculator(m, 2);

        for (LPOrder l = 0; l < m_nstr; ++l) {
            for (StreamIndex i = 0; i < m_nstr; ++i) {
                double theta = acos(m_abscissae[i]);

                m_lp_abscissae4[m][i][l].P() = calculatorP.d(theta, l);
                m_lp_abscissae1[m][i][l].value = calculatorP.d(theta, l);

                m_lp_abscissae4[m][i][l].R() =
                    -0.5 *
                    (calculatorpos.d(theta, l) + calculatorneg.d(theta, l));
                m_lp_abscissae4[m][i][l].T() =
                    -0.5 *
                    (calculatorpos.d(theta, l) - calculatorneg.d(theta, l));

                if (m % 2 == 1) {
                    // m_lp_abscissae4[m][i][l].R() *= -1;
                    // m_lp_abscissae4[m][i][l].T() *= -1;
                }
            }
        }
    }
}

sasktran_disco::SKTRAN_DO_UserSpec*
sasktran_disco::SKTRAN_DO_UserSpec::setNumberOfLayers(
    unsigned int num_atmo_layers) {
    using namespace sasktran_disco;
    if (num_atmo_layers < 1) {
        throw InvalidConfiguration(
            "Number of layers must be greater than zero!");
    } else {
        m_nlyr = num_atmo_layers;
    }
    return this;
}

sasktran_disco::SKTRAN_DO_UserSpec*
sasktran_disco::SKTRAN_DO_UserSpec::setTOAIntensities(double direct) {
    using namespace sasktran_disco;
    if (direct < 0) {
        throw InvalidConfiguration(
            "Radiances at the top of the atmosphere cannot be negative!");
    } else {
        m_itop_direct = direct;
    }
    return this;
}

sasktran_disco::SKTRAN_DO_UserSpec*
sasktran_disco::SKTRAN_DO_UserSpec::setSSAEqual1Dither(double dither) {
    if (dither <= 0) {
        using namespace sasktran_disco;
        throw InvalidConfiguration("SSA is 1 dither must greater than zero!");
    }
    m_ssalb1_dither = dither;
    return this;
}

#pragma endregion

#pragma region "Getters"
sasktran_disco::uint
sasktran_disco::SKTRAN_DO_UserSpec::getNumberOfStreams() const {
    return m_nstr;
}

sasktran_disco::uint
sasktran_disco::SKTRAN_DO_UserSpec::getNumberOfLayers() const {
    return m_nlyr;
}

double sasktran_disco::SKTRAN_DO_UserSpec::getTopDirectIntensity() const {
    return m_itop_direct;
}

const sasktran_disco::VectorDim1<double>*
sasktran_disco::SKTRAN_DO_UserSpec::getStreamAbscissae() const {
    return &m_abscissae;
}
const sasktran_disco::VectorDim1<double>*
sasktran_disco::SKTRAN_DO_UserSpec::getStreamWeights() const {
    return &m_weights;
}
const sasktran_disco::VectorDim3<sasktran_disco::LegendrePhaseContainer<4>>*
sasktran_disco::SKTRAN_DO_UserSpec::getAbscissaeLegendreP4() const {
    return &m_lp_abscissae4;
}

const sasktran_disco::VectorDim3<sasktran_disco::LegendrePhaseContainer<1>>*
sasktran_disco::SKTRAN_DO_UserSpec::getAbscissaeLegendreP1() const {
    return &m_lp_abscissae1;
}

double sasktran_disco::SKTRAN_DO_UserSpec::getSSAEqual1Dither() const {
    return m_ssalb1_dither;
}

#pragma endregion
