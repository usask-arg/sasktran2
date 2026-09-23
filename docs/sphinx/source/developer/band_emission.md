(_dev_band_emission)=
# Extending band emission

Band emission separates an independent photon VER profile from the model that
distributes those photons among spectral lines. The native extension point is
`sasktran2_rs::emission::BandEmissionModel`. The Python `O2BandEmissionRate`
interface selects the O2 implementation of this framework.

There are three parts:

- `BandLineData` contains the band name, line wavelengths in nm, and molecular
  mass in atomic mass units. Its constructor validates the wavelengths and mass.
- `BandEmissionModel` owns the spectroscopy and population assumptions used to
  calculate normalized relative line intensities. Einstein coefficients, level
  energies, and independent excitation parameters belong here.
- `BandVolumeEmissionRate<M>` owns the altitude and photon VER profiles. It
  handles interpolation, isotropic emission, Doppler broadening, spectral-grid
  conversion, and derivative mappings using model `M`.

The generic source has no molecular selection logic, Einstein coefficients, or
assumption that rotational populations follow the kinetic temperature. The O2
implementation in `emission/o2.rs` supplies its molecular mass, selects the
vibrational transition, and evaluates the existing rotational LTE weight models.
Population-to-VER conversion remains in `PopulationEmissionRate`.

## Implementing a population model

Implement these methods in Rust:

```rust
pub trait BandEmissionModel {
    fn line_data(&self) -> &BandLineData;

    fn line_weights(&self, inputs: &impl StorageInputs) -> Result<Array2<f64>>;

    fn line_weights_with_temperature_derivative(
        &self,
        inputs: &impl StorageInputs,
    ) -> Result<(Array2<f64>, Array2<f64>)>;
}
```

Weights have shape `(atmospheric location, line)`, with one column per wavelength
in `BandLineData`, in the same order. Each row must be finite, nonnegative, and
sum to one. Locations are those exposed by `StorageInputs`; a 2D atmosphere uses
its flattened native locations. Independent model profiles must be evaluated or
interpolated onto these locations before returning weights. The common source
checks array shapes; each model is responsible for validating and normalizing
its weights.

The derivative method returns the weights and their **local partial derivative
with respect to atmospheric kinetic temperature**, at fixed total band VER and
fixed independent excitation parameters. Derivative rows must sum to zero.
These derivatives exclude Doppler broadening, which the source calculates using
the kinetic temperature and the supplied mass. Temperature-independent weights
therefore return zero weight derivatives, while the emitted spectrum still has
a Doppler temperature derivative. A population model must explicitly provide
this method; missing derivatives are not silently treated as zero.

Only the value method is called on the path without temperature derivatives.
Keep derivative allocation and evaluation out of that method. The source uses
static dispatch, so the native O2 path introduces no Python callbacks or virtual
dispatch in the line-profile loop.

Construct `BandVolumeEmissionRate::new(altitudes_m, photon_ver, model)` with the
new model. VER is in photons m^-3 s^-1 integrated over all band lines and solid
angle. The existing VER and atmospheric-temperature derivative mappings are
shared by all implementations. Separate retrieval parameters for excitation
temperatures or individual level populations require additional mappings; they
are not exposed by this interface yet. Emission profiles currently use Doppler
broadening. Absorption remains a separate constituent.

## Validation

The generic source tests use synthetic line populations, independent of O2
spectroscopy. They check photon conservation, mass-dependent widths, analytic
temperature derivatives for both fixed and temperature-dependent populations,
and that the value-only path never requests population derivatives. O2 tests
add radiance-level VER and temperature finite differences, including
self-absorption, population conversion, alternate altitude grids, and both
line-weight conventions.
