---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
  language: python
---

(_example_stratospheric_aerosol)=
# Stratospheric aerosol reference scenarios

`sk.climatology.stratospheric_aerosol` provides twelve fixed sulfate aerosol
reference cases derived from the USask SAGE III–ISS particle-size retrieval.
Each case keeps an observed extinction profile paired with its retrieved median
radius. The small catalogue is bundled with SASKTRAN2; the original archive is
not needed to select or prepare profiles.

```{code-cell}
import numpy as np
import matplotlib.pyplot as plt
import sasktran2 as sk

aerosol = sk.climatology.stratospheric_aerosol
aerosol.scenarios()[[
    "event_id", "latitude_band", "loading", "loading_percentile",
    "observed_bottom_m", "observed_top_m", "reference_upper_scale_height_m",
]]
```

The latitude bands are southern midlatitudes (55–35°S, `sh_midlat`), tropics
(20°S–20°N, `tropical`), and northern midlatitudes (35–55°N, `nh_midlat`). Each
has `low`, `typical`, `elevated`, and `extreme` cases, selected near the 10th,
50th, 90th, and 99th percentiles of **756 nm optical depth over 18–30 km**.
Percentiles use the screened June 2017–July 2026 archive within each band.
They describe relative loading in that sample, not global severity classes.
The two polar ERS latitude bands have insufficient SAGE sampling for this catalogue.

Selection uses actual paired observations close to the median extinction/radius
shape within each percentile neighbourhood. Source aerosol flags 2/3 are required
in all four near-infrared channels, alongside positive extinction, bounded radius,
formal relative errors below 50%, a tropopause margin, and spectral consistency.
The actual contiguous reliable interval is retained, including valid measurements
below 18 km or above 30 km. Internal gaps are not filled.

```{code-cell}
fig, axes = plt.subplots(1, 3, figsize=(12, 5), sharey=True)
bands = {"sh_midlat": "Southern midlatitudes", "tropical": "Tropics",
         "nh_midlat": "Northern midlatitudes"}
for ax, (band, title) in zip(axes, bands.items()):
    for tier in ["low", "typical", "elevated", "extreme"]:
        p = aerosol.profile(f"{band}_{tier}")
        core = p.region == 0
        line, = ax.semilogx(
            p.extinction_per_m.where(core), p.altitude_m / 1000, label=tier
        )
        ax.semilogx(
            p.extinction_per_m.where(
                ~core | (p.altitude_m == p.observed_bottom_m)
                | (p.altitude_m == p.observed_top_m)
            ), p.altitude_m / 1000,
            "--", color=line.get_color(), alpha=0.65,
        )
    ax.set(title=title, xlabel="Extinction at 756 nm [m$^{-1}$]",
           xlim=(1e-11, 1e-4), ylim=(0, 55))
    ax.grid(alpha=0.2)
axes[0].set_ylabel("Altitude [km]")
axes[-1].legend()
fig.tight_layout()
```

Solid curves show the observed interval after smoothing. Dashed curves show
modelled extensions. The defaults are:

- Gaussian smoothing of **log extinction**, with 1.5 km full width at half maximum,
  followed by rescaling to conserve optical depth over the observed interval.
  Radius is already regularized in the source retrieval and is left unchanged.
- A downward exponential taper with scale height 2 km, modified to reach exactly
  zero at the ground. This represents the lower tail of this stratospheric
  constituent; add a separate constituent if tropospheric aerosol is needed.
- An upper exponential continuation with a fixed regular-condition extinction
  scale height for each latitude band. All four tiers share that rate; their
  extinction at the joining altitude remains their own. A different smoothing
  width or loading tier does not refit the upper scale height.
- Constant endpoint median radius in each extension, with lognormal width 1.6.

The upper reference is calibrated from the 10th–60th loading percentiles,
requiring background flag 2 throughout 26–30 km. Robust raw log-extinction slopes
over 27–30 km are aggregated by month and then across months, with at least five
profiles per contributing month. The scale heights are rounded to 100 m; the
bundled build report records alternative fit intervals and smoothing sensitivity.
Raw fitting avoids a reflected smoothing boundary altering the calibration.

| Latitude band | Default upper scale height | Raw-fit sensitivity across the three tested intervals |
| --- | ---: | ---: |
| Southern midlatitudes | 2.8 km | 2.83–3.00 km |
| Tropics | 3.6 km | 3.55–4.42 km |
| Northern midlatitudes | 2.8 km | 2.76–3.47 km |

The sensitivity ranges describe method dependence, not confidence intervals.

Published comparisons include a mean 3.2 km extinction scale height in background
SAGE profiles ([Brogniez and Lenoble, 1987](https://doi.org/10.1029/JD092iD03p03051)),
3.75 km above 26 km during volcanic aerosol abatement
([Elterman et al., 1969](https://doi.org/10.1364/AO.8.000893)), and a 4 km
continuation above 30 km in a SCIAMACHY retrieval discussion paper
([Ernst et al., 2012, §3.4](https://amt.copernicus.org/preprints/5/5993/2012/amtd-5-5993-2012-print.pdf)).
These are comparisons, not universal bounds or an uncertainty interval. Neither
these references nor our fitted rates validate constant sulfate size and scale
height all the way to 100 km; the highest levels are a numerical continuation.

## Inspect and customize a profile

```{code-cell}
p = aerosol.profile("tropical_extreme")
p[["core_aod", "lower_extension_aod", "upper_extension_aod_to_infinity"]]
```

Those three optical depths separate the observed-core contribution from the two
modelled additions. They describe the continuous prepared profile independently
of the output grid. The upper integral extends to infinity; the default grid ends
at 100 km. Resampling to a coarse grid can change the numerical integral on that
grid. Loading labels are fixed before adding extensions.

```{code-cell}
custom = aerosol.profile(
    "tropical_extreme",
    altitudes_m=np.arange(0., 65001., 250.),
    smoothing_fwhm_m=1000.,
    lower_scale_height_m=1500.,
    upper_scale_height_m=3200.,  # 4000 is another literature comparison
)

unprocessed = aerosol.profile(
    "tropical_extreme",
    smoothing_fwhm_m=0.,
    lower_extension="zero",
    upper_extension="zero",
)

raw = aerosol.load_dataset()
```

`load_dataset()` retains all nine measured extinction channels and their formal
errors, median radius and its formal error, source flags, actual wavelengths,
event identity and source-file checksums. `observed_valid` identifies the usable
interval of each case on the shared native altitude coordinate. Padding outside
that interval is missing data. Prepared profiles retain raw variables on
`observed_altitude_m`; model outputs use `altitude_m` and a `region` flag (-1
lower extension, 0 observed core, 1 upper extension). The original formal errors
are not propagated through smoothing and are not uncertainties of the tails.

For an offline or modified catalogue, all helpers accept `path="catalogue.nc"`.
This bypasses the pinned checksum and records the actual file checksum. Its schema,
units and selected profile are still validated. Alternatively use `db_root` to
choose where the bundled catalogue is cached. The two options are mutually exclusive.

## Add to an atmosphere

The aerosol helper is independent of ERS gas, temperature and pressure selection:

```python
atmosphere["stratospheric_aerosol"] = aerosol.constituent(
    "tropical_typical",
    altitudes_m=atmosphere.model_geometry.altitudes(),
)
```

This returns the existing `ExtinctionScatterer`, normalized at 756 nm, with sulfate
Mie scattering, the paired median radius, and fixed width 1.6. Passing the model
altitudes evaluates the extensions directly on that grid. Otherwise the default
prepared grid is used; the constituent is zero outside its grid. The default Mie
property calculates optics on demand and may download the standard OSIRIS H2SO4
refractive-index file on first use.

For repeated or large spectral calculations, supply a compatible cached table:

```python
optics = sk.database.MieDatabase(
    sk.mie.LogNormalDistribution().freeze(mode_width=1.6),
    sk.mie.refractive.H2SO4(),
    np.unique(np.r_[atmosphere.wavelengths_nm, 756.]),
    median_radius=np.arange(10., 600., 10.),
)
atmosphere["stratospheric_aerosol"] = aerosol.constituent(
    "tropical_typical", optics,
    altitudes_m=atmosphere.model_geometry.altitudes(),
)
```

The caller is responsible for the composition, fixed width, wavelength and radius
coverage of custom optics. The fixed-width sulfate retrieval does not establish
smoke absorption. These profiles also do not imply a season or local time matching
an independently chosen gas scenario.
