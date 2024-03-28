#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_lazyazimuth.h"
#include <sasktran2/math/wigner.h>

void sasktran_disco::AlbedoExpansion::calculateAEOrder(
    AEOrder m, sasktran_disco::Albedo& albedo) {
    // Configure the BRDF object for the current order of the azimuth expansion.
    albedo.configure(m, M_LOS, M_MU, M_CSZ, m_brdf.get(), m_nterms);
}
