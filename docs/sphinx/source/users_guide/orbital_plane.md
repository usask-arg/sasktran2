# Orbital-plane 2D calculations

`OrbitalPlaneEngine` evaluates timestamped ECEF limb measurements against a
full along-track atmosphere. It divides measurements into fixed, half-open time
groups, fits a local great-circle plane to each group, retains every native 2D
engine, and stitches radiance and derivative results back into the original
line-of-sight order.

The atmospheric layout is `(orbital_position, altitude)`, flattened with
altitude varying fastest internally. Pressure, temperature, specific humidity,
and refractive index may be supplied either as one altitude profile or on the
full orbital grid. The 2D absorber and scatterer constituents use the native
orbital layout directly.

```python
viewing = sk.OrbitalPlaneViewingGeometry.from_tangent_locations(
    times,
    observer_positions_ecef_m,
    tangent_locations_ecef_m,
    vertical_slice=vertical_slice,
    geoid=sk.WGS84(),
)
geometry = viewing.construct_atmosphere_geometry(
    altitude_grid_m,
    along_track_angle_delta=0.01,  # radians
    path_padding_angle=np.deg2rad(10),
    max_orbital_positions=100_000,
)

config.los_refraction = True
engine = sk.OrbitalPlaneEngine(
    config,
    geometry,
    viewing,
    time_group_duration_s=30,
    group_padding_angle=np.deg2rad(10),
    max_horizontal_scale_residual=0.01,
    solar_handler=solar_handler,
    derivative_execution="resident",
    max_eager_jacobian_bytes=2 * 1024**3,
)

atmosphere = sk.Atmosphere(
    geometry,
    config,
    wavelengths_nm=wavelengths_nm,
    calculate_derivatives=False,
)
atmosphere.pressure_pa = pressure_pa
atmosphere.temperature_k = temperature_k
result = engine.calculate_radiance(atmosphere)
```

A derivative-enabled atmosphere does not force every forward call to allocate
weighting functions. Use
``engine.calculate_radiance(atmosphere, derivatives=False)`` to obtain the same
radiance (and requested LOS optical depth) without group derivative workspaces
or a stitched Jacobian. The atmosphere remains derivative enabled and can be
passed to ``linearize()`` afterward.

By default, both forward and derivative calculations use resident group
engines. This is intended for retrievals: after the first calculation, optical
state updates reuse the local engines, atmosphere allocations, traced paths,
and integration workspaces. Set ``derivative_execution="streaming"`` to create
one temporary derivative engine at a time during a VJP when lower retained
memory is more important than repeated-call speed. Forward calculations remain
resident in either mode.

Observer positions and look directions must have shape ``(los, 3)`` and
timestamps must have shape ``(los,)`` without ``NaT``. ``vertical_slice`` is a
required integer identifier with shape ``(los,)``. All LOS from one limb image
or vertical scan use the same identifier and occupy one contiguous block. This
flat representation supports both fixed-size images and ragged vertical scans.
Look directions are normalized once on construction; the normalized values
exposed by the viewing geometry are read-only and are the same values used
natively. Input arrays are copied, so construction never normalizes a
caller-owned array in place.

Track construction uses one spherical-mean tangent-location anchor per vertical
slice. Individual LOS order within an image or scan therefore does not affect
the master atmosphere grid. Every LOS still retains its own time, observer,
look direction, tangent altitude, and exact transformed horizontal location.
For a finite-altitude limb image, a least-squares tangent-location gradient with
tangent altitude supplies the local track direction independently of row order.
A timestamp gradient supplies the direction for a true sequential vertical scan;
only a geometrically degenerate image needs the boresight fallback. The full
angular extent of each slice is included before path padding is added.

Time grouping is also slice atomic. Each complete slice is assigned to the
half-open bin containing its mean sample time, while the group's reference time
is the mean of all actual LOS samples assigned to that group.

For an operational retrieval, the orbit can also be divided into a small
number of contiguous vertical-slice chunks (for example, four), with an
``OrbitalPlaneEngine`` retained for each chunk. Radiances concatenate in the
original observation order and VJP contributions scatter-add onto the shared
master atmospheric grid. Chunking is not required by the interface; it is an
outer memory/throughput choice.

```python
viewing_chunks = viewing.split(4, time_group_duration_s=30)
chunk_engines = [
    sk.OrbitalPlaneEngine(
        config,
        geometry,  # one shared master atmosphere grid
        chunk,
        time_group_duration_s=30,
        solar_handler=solar_handler,
    )
    for chunk in viewing_chunks
]
```

``viewing.split(...)`` requires the same group duration used by the engines. It
never divides an image, vertical scan, or native time group, and every chunk
retains the original full-orbit time-bin origin. ``viewing.isel`` is also
available for positional LOS subsets, provided each retained slice identifier
remains one contiguous block.

The constructed master grid follows the selected geoid at zero altitude and
preserves each sampled normalized ECEF direction by intersecting that radial
line directly with the reference ellipsoid. It is
extended beyond the first and last nominal tangent locations by
`path_padding_angle` on each end. The default is 10 degrees. If a boundary
group's paths need more atmosphere than this,
`edge_clipping` in `group_diagnostics` reports the clamped edge.
`geometry.surface_radii_m` contains the geocentric geoid radius at every
orbital position. Each time group's native 2D calculation remains locally
spherical. Its constant radius is the mean interpolated geoid radius at the
actual observation tangent locations, which is the least-squares choice for
the local conversion from angular separation to horizontal distance. The value
is reported as `earth_radius_m` in `group_diagnostics`. The wider padding
window does not affect this radius.

The original ECEF LOS is not passed unchanged into that approximate sphere.
Rust first derives the LOS's tangent altitude, tangent latitude/longitude, and
observer altitude relative to the selected geoid. It then expresses the geoid
surface location in the fitted group's horizontal and out-of-plane angles and
reconstructs the ray with the geometry-relative tangent-altitude policy. Thus
the vertical tangent coordinate and horizontal geoid location are conserved
explicitly, while the single group radius only controls the residual path and
horizontal-distance approximation. The per-ray policy values and the remaining
surface-radius scale residual are available in `group_diagnostics`.
Independently, every retained group expands the interval between its first and
last tangent locations by ``group_padding_angle`` on each side. No additional
horizon margin is added. This margin is user configurable and defaults to 10
degrees. ``padding_angle`` in the group diagnostics reports the applied value.
At a physical end of the master grid,
the window still clamps and reports ``edge_clipping``; increase the master
``path_padding_angle`` as well when those boundary groups need the full margin.
Atmospheric altitudes are interpreted as heights above the geoid. For a
spherical model, the explicit `OrbitalPlaneGeometry(earth_radius_m, ...)`
constructor remains available.
Zero path padding is supported. For a single unique tangent location, the
constructor still creates the minimum three-point track needed to define a 2D
plane. ``max_orbital_positions`` protects against accidentally requesting an
enormous grid with a very small angular spacing; set it to ``None`` only when a
larger allocation is intentional.

The constructed grid carries its mapping back to the measurements. The
``geometry.grid_dataset`` dataset contains the ECEF surface point, surface
radius, cumulative along-track angle, geodetic latitude/longitude, and an
interpolated reference time for every ``orbital_position``. Reference time is
linearly interpolated between vertical-slice mean times and clamped through the
two padding regions. ``geometry.vertical_slice_anchors`` retains each original
slice identifier, mean timestamp, ECEF anchor, and along-track coordinate.
These auxiliary coordinates are also attached to orbital-position retrieval
parameters in ``linearization.tangent_template`` and the eager Jacobian.

Ancillary locations can be mapped onto the same grid without reimplementing
the track logic:

```python
mapping = geometry.project_ground_track(
    ancillary_ground_locations_ecef_m,
    order="increasing",
)
ancillary_on_grid = np.interp(
    geometry.cumulative_angles,
    mapping.along_track_angle_rad,
    ancillary_values,
)
```

The result includes the containing segment, its interpolation fraction, the
along-track coordinate, cross-track angular residual, and an edge flag.
``order="increasing"`` or ``"decreasing"`` constrains successive searches to
the remaining ordered track, disambiguating closed-orbit endpoints and
self-crossings. Use ``"independent"`` only when the inputs have no meaningful
sequence.

When line-of-sight refraction is enabled, an explicit
`atmosphere.refractive_index` takes precedence. Otherwise the engine evaluates
Ciddor refractivity in Rust at 600 nm and 400 ppm CO2 by default, treating
missing humidity as zero. Each ray interpolates a locally spherical vertical
profile with the same projected local-2D horizontal cell and fraction used for
its optical atmosphere. Exact solar paths remain straight, and
`solar_refraction=True` is rejected.

Use `engine.group_diagnostics` to inspect each group's observation and master
grid indices, reference time and basis, local mean Earth radius, projected
horizontal angles and distances, tangent-altitude and horizontal-scale
residuals, plane-fit residual, edge clipping, traced-path edge usage, window
expansion status, geometry refresh count, and persistent engine/workspace
identities.
``maximum_relative_horizontal_scale_residual`` measures the complete local
horizontal-distance error, including fitted-plane angular distortion and the
group's mean-radius approximation over the entire selected window. The default
``max_horizontal_scale_residual=0.01`` rejects a group above 1%; set it to
``None`` only when accepting that approximation explicitly. Separate angular,
radius, combined-distance, and observation-radius diagnostics are also exposed.

Dynamic post-refraction window expansion is not needed for the normal path:
the fixed ``group_padding_angle`` is included when each internal engine is
constructed and remains resident with it. ``window_expanded`` therefore stays
false in this version. A path that exceeds even that margin uses the structured
grid's extended horizontal edge cell. ``edge_clipping`` reports truncation of
the requested group window at the master grid, while ``path_edge_clipping`` and
``path_edge_clipping_per_los`` report actual traced atmospheric segments beyond
the local grid, including after refraction changes. Choose the group and master
padding large enough to keep both kinds of clipping false for retrieval work.

Configuration values that determine native engine structure are frozen at
engine construction. Mutating sources, Stokes dimension, refraction mode,
threading structure, or related settings on the original ``Config`` raises a
clear error on the next calculation; construct a new orbital engine after such
a change. Atmospheric optical and thermodynamic state remains mutable and is
the intended retrieval state.

`linearize()` freezes the current ray paths for its lazy Jacobian, JVP, and VJP.
Pressure, temperature, and humidity derivatives therefore contain optical-state
effects but do not differentiate motion of a refracted path. Altitude-only
inputs have an `altitude` parameter dimension; native fields have
`(orbital_position, altitude)` dimensions. Explicit refractive index is
geometry-only in this version. Installing different refractive geometry or a
different spatial Lambertian field in the same engine invalidates any older
linearization, which then raises ``StaleLinearizationError`` instead of silently
using the newer state.

An eager full orbital Jacobian can be much larger than a VJP. By default,
``calculate_radiance`` with a derivative-enabled atmosphere and the lazy
``linearization.jacobian`` property are rejected when their stitched derivative
payload is estimated to exceed 2 GiB. Use
``calculate_radiance(..., derivatives=False)`` for radiance-only work or
``linearize()`` followed by JVP/VJP for retrievals. Advanced callers can inspect
``engine.estimate_eager_jacobian_bytes(atmosphere)`` and explicitly disable the
guard with ``max_eager_jacobian_bytes=None``.

For JVP and VJP products, Rust copies only the derivative mappings named by the
supplied tangent or requested VJP parameters into each local group atmosphere.
Repeated calls with the same selection reuse that resident allocation. The
active names are exposed as ``resident_volume_derivative_mappings`` and
``resident_surface_derivative_mappings`` in ``group_diagnostics``. This makes a
one-parameter or small-block retrieval independent of unrelated mappings that
remain registered on the master atmosphere.

``NumberDensityScatterer2D`` and ``ExtinctionScatterer2D`` preserve the input
topology of optical-property arguments such as aerosol radius: scalar inputs
produce scalar retrieval parameters, ``(altitude,)`` inputs produce altitude
profiles shared along track, and ``(orbital_position, altitude)`` inputs remain
fully two dimensional. Optical calculations broadcast the first two forms at
the native boundary; VJPs sum the local contributions back to the supplied
parameter grid.

The orbital engine supports exact single scattering, successive-orders
multiple scattering, occultation/transmission, standard emission,
volume-emission-rate calculations, limb rays, and surface-intersecting ECEF
rays. Flux observers, diffuse-ray refraction, and solar-path refraction are not
supported. `LambertianSurface2D` accepts scalar,
gray along-track, and along-track/spectral albedo fields. Rust gathers the
actual group rows without averaging, and the native ground source linearly
interpolates albedo at each traced surface intersection. Its eager Jacobian,
JVP, and VJP use the same current traced stencil, including after a refractive
geometry refresh. In successive-orders calculations, discrete ground source
points use a unit Lambertian albedo; every ground-terminating transport ray
then applies the interpolated albedo at its physical surface intersection.
