#pragma once

#include "atmosphere.h"
#include "sasktran2/atmosphere/surface.h"
#include "sasktran2/derivative_mapping.h"
#include "sasktran2/output.h"
#include <memory>
#include <sasktran2.h>

struct Config {
    sasktran2::Config impl;

    Config() {}
};

struct Geodetic {
    std::unique_ptr<sasktran2::math::geodetic::Geodetic> impl;

    Geodetic(double equatorial_radius, double flattening_factor)
        : impl(std::make_unique<sasktran2::math::geodetic::Geodetic>(
              equatorial_radius, flattening_factor)){};
};

struct DerivativeMapping {
    sasktran2::DerivativeMapping* impl;

    DerivativeMapping(sasktran2::DerivativeMapping* mapping) : impl(mapping){};
};

struct SurfaceDerivativeMapping {
    sasktran2::SurfaceDerivativeMapping* impl;

    SurfaceDerivativeMapping(sasktran2::SurfaceDerivativeMapping* mapping)
        : impl(mapping){};
};

struct BRDF {
    std::shared_ptr<sasktran2::atmosphere::brdf::BRDFInterface> impl;
};

struct Geometry1D {
    std::unique_ptr<sasktran2::Geometry1D> impl;

    Geometry1D(double cos_sza, double saa, double earth_radius,
               double* grid_values, int ngrid_values, int interp_method,
               int geotype);
};

struct AtmosphereStorage {
    std::unique_ptr<sasktran2::atmosphere::AtmosphereGridStorage> impl;

    AtmosphereStorage(int nlocation, int nwavel, int nphase_moments,
                      int nstokes, double* ssa, double* total_extinction,
                      double* emission_source, double* leg_coeff,
                      double* solar_irradiance);

    int get_derivative_mapping(const char* name, DerivativeMapping** mapping);
};

struct Surface {
    std::unique_ptr<sasktran2::atmosphere::SurfaceInterface> impl;

    Surface(int nwavel, int nstokes, double* emission);
};

struct Atmosphere {
    std::unique_ptr<sasktran2::atmosphere::AtmosphereInterface> impl;

    Atmosphere(AtmosphereStorage* storage, Surface* surface,
               bool calculate_derivatives, bool calculate_emission_derivatives);
};

struct OutputC {
    std::unique_ptr<sasktran2::OutputInterface> impl;

    OutputC(double* radiance, int nrad, int nstokes);

    int assign_derivative_memory(const char* name, double* derivative_mapping,
                                 int nrad, int nstokes, int nderiv);
    int assign_surface_derivative_memory(const char* name,
                                         double* derivative_mapping, int nrad,
                                         int nstokes);
};

struct ViewingGeometry {
    sasktran2::viewinggeometry::ViewingGeometryContainer impl;

    ViewingGeometry() {}
};
