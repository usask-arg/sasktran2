#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_pconfig.h"
#include "sktran_disco/sktran_do_polarization_types.h"
#include "sktran_disco/sktran_do_specs.h"
#include "sktran_disco/sktran_do_types.h"

template <int NSTOKES, int CNSTR>
void sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR>::configure(
    SKTRAN_DO_UserSpec& userspecmemory, const sasktran2::Config& config,
    double cos_sza, int nlyr,
    const std::vector<sasktran2::raytracing::TracedRay>& traced_rays) {

    int nstr = config.num_do_streams();

    m_userspec = &userspecmemory;

    const_cast<double&>(this->M_CSZ) = cos_sza;
    const_cast<double&>(this->M_SAZ) = 0;
    const_cast<bool&>(this->M_SOLAR_REFRACTION) = config.solar_refraction();

    const_cast<uint&>(this->M_NSTR) = config.num_do_streams();
    const_cast<uint&>(this->M_NLYR) = nlyr;

    userspecmemory.configure(nstr, nlyr);

    this->M_MU = m_userspec->getStreamAbscissae();
    this->M_WT = m_userspec->getStreamWeights();
    configureLP(m_userspec);

    const_cast<bool&>(this->M_SS_ONLY) = false;
    const_cast<bool&>(this->M_BACKPROP_BVP) = config.do_backprop();

    m_lp_csz_storage.resize(
        this->M_NSTR,
        VectorDim1<LegendrePhaseContainer<NSTOKES>>(this->M_NSTR));

    this->M_LP_CSZ = &m_lp_csz_storage;
    m_pool.init(this->M_NLYR, this->M_NSTR);
    m_poolindex = 0;

    for (int m = 0; m < (int)this->M_NSTR; ++m) {
        for (int l = 0; l < (int)this->M_NSTR; ++l) {
            m_lp_csz_storage[m][l].fill(m, l, this->M_CSZ);
        }
    }
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::PersistentConfiguration<
    NSTOKES, CNSTR>::configureModelSpecs(const SKTRAN_DO_UserSpec* userspec) {
    m_userspec = userspec;
    const_cast<uint&>(this->M_NSTR) = m_userspec->getNumberOfStreams();
    const_cast<uint&>(this->M_NLYR) = m_userspec->getNumberOfLayers();
    this->M_MU = m_userspec->getStreamAbscissae();
    this->M_WT = m_userspec->getStreamWeights();
    configureLP(userspec);
    const_cast<bool&>(this->M_SS_ONLY) = m_userspec->getSSOnly();
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR>::configureLP(
    const SKTRAN_DO_UserSpec* userspec) {
    if constexpr (NSTOKES == 1) {
        this->M_LP_MU = userspec->getAbscissaeLegendreP1();
    } else {
        this->M_LP_MU = (const VectorDim3<LegendrePhaseContainer<NSTOKES>>*)
                            userspec->getAbscissaeLegendreP4();
    }
}

SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(sasktran_disco::PersistentConfiguration);
