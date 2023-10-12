import sasktran as sk
import numpy as np
from sasktran.mie.distribution import sp_lognormal, integrated_mie
from sasktran.mie.refractive import refractive_index_fn_h2so4, refractive_index_fn_water
from sasktran.legendre import compute_greek_coefficients_legendre


def create_table(name: str, distribution, refractive_index_fn, wavelengths):
    mie = sk.MieWiscombe()

    vals = integrated_mie(mie, distribution, refractive_index_fn, wavelengths, num_quad=1000, maxintquantile=0.99999)

    # vals xs is in units of nm^2, convert to cm^2
    vals['xs_total'] *= 1e-14
    vals['xs_scattering'] *= 1e-14
    vals['xs_absorption'] *= 1e-14

    lm_a1 = np.zeros_like(vals['lm_p11'].values)
    lm_a2 = np.zeros_like(lm_a1)
    lm_a3 = np.zeros_like(lm_a1)
    lm_a4 = np.zeros_like(lm_a1)
    lm_b1 = np.zeros_like(lm_a1)
    lm_b2 = np.zeros_like(lm_a1)

    for idx in range(len(vals.wavelength.values)):
        selected = vals.isel(wavelength=idx)

        lm_a1[idx, :], lm_a2[idx, :], lm_a3[idx, :], lm_a4[idx, :], lm_b1[idx, :], lm_b2[idx, :] = compute_greek_coefficients_legendre(
            selected['lm_p11'].values, selected['lm_p12'].values, selected['lm_p11'].values,
            selected['lm_p33'].values, selected['lm_p34'].values, selected['lm_p33'].values
        )

    vals['lm_a1'] = (['wavelength', 'legendre'], lm_a1)
    vals['lm_a2'] = (['wavelength', 'legendre'], lm_a2)
    vals['lm_a3'] = (['wavelength', 'legendre'], lm_a3)
    vals['lm_a4'] = (['wavelength', 'legendre'], lm_a4)
    vals['lm_b1'] = (['wavelength', 'legendre'], lm_b1)
    vals['lm_b2'] = (['wavelength', 'legendre'], lm_b2)

    return vals

database_wavelength = np.arange(250, 1500, 5.0)

# Create the table for every entry
# SULPHATE - FINE AEROSOL
dist = sp_lognormal(40, 1.3)
refrac = refractive_index_fn_h2so4()

vals = create_table('sulfate_fine', dist, refrac, database_wavelength)

# SULPHATE - COARSE AEROSOL
dist = sp_lognormal(120, 1.3)
refrac = refractive_index_fn_h2so4()

create_table('sulfate_coarse', dist, refrac, database_wavelength)

# SULPHATE - STRATOSPHERIC AEROSOL
dist = sp_lognormal(80, 1.6)
refrac = refractive_index_fn_h2so4()

create_table('sulfate_strat', dist, refrac, database_wavelength)

# WATER CLOUD
dist = sp_lognormal(8000 / np.exp(5/2 * np.log(1.2)**2), 1.2)
refrac = refractive_index_fn_water()
create_table('water_cloud', dist, refrac, database_wavelength)
