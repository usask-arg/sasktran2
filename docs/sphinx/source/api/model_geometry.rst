.. _api_model_geometry:

Model Geometry
==============

.. autosummary::
    :toctree: generated/

    sasktran2.Geometry1D
    sasktran2.Geometry2D

Geometry2D capabilities
-----------------------

The current :class:`sasktran2.Engine` integration for
:class:`sasktran2.Geometry2D` supports occultation and transmission, exact
single scattering, standard emission, and volume-emission-rate sources.
Altitude-only constituents remain compatible and are applied across the
horizontal grid, while the native 2D constituents can specify horizontally
varying volume properties. Existing surface constituents may be used, but the
surface itself is spatially uniform.

Multiple scattering, flux observers, line-of-sight refraction, and
horizontally varying surfaces are not yet supported with
:class:`sasktran2.Geometry2D`.

Geometry2D engine calculations use the Rust-backed 2D ray tracer and therefore
require Rust component support. The ``SKTRAN_USE_RUST_RAYTRACER`` CMake option
selects the Rust-backed ray tracer only for Geometry1D; disabling it does not
disable the Geometry2D ray tracer.
