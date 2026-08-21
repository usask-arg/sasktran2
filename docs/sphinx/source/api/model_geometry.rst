.. _api_model_geometry:

Model Geometry
==============

.. autosummary::
    :toctree: generated/

    sasktran2.Geometry1D
    sasktran2.Geometry2D
    sasktran2.OrbitalPlaneGeometry
    sasktran2.OrbitalPlaneViewingGeometry

Geometry2D capabilities
-----------------------

The current :class:`sasktran2.Engine` integration for
:class:`sasktran2.Geometry2D` supports occultation and transmission, exact and
table single scattering, successive-orders multiple scattering, standard
emission, and volume-emission-rate sources. Altitude-only constituents remain
compatible and are applied across the horizontal grid, while native 2D
constituents can specify horizontally varying volume properties. Existing
surface constituents may be used, but the surface itself is spatially uniform.

Solar-ray refraction is supported by the exact and table single-scatter
sources and by the direct-beam calculation used by successive orders. Other
multiple-scatter sources, flux observers, successive-orders diffuse-ray
refraction, and a public horizontally varying surface state are not supported
by the ordinary :class:`sasktran2.Geometry2D` engine.
The :class:`sasktran2.OrbitalPlaneEngine` adds per-ray line-of-sight refraction
and a fully along-track/spectral Lambertian field, interpolated at each traced
surface intersection, while retaining local Geometry2D engines.

Geometry2D engine calculations use the Rust-backed 2D ray tracer and therefore
require Rust component support. The ``SKTRAN_USE_RUST_RAYTRACER`` CMake option
selects the Rust-backed ray tracer only for Geometry1D; disabling it does not
disable the Geometry2D ray tracer.
