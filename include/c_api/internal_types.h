#pragma once

#include "sasktran2/derivative_mapping.h"
#include "sasktran2/output.h"
#include <sasktran2.h>

struct Config {
    sasktran2::Config impl;

    Config() {}
};

struct DerivativeMapping {
    sasktran2::DerivativeMapping* impl;

    DerivativeMapping(sasktran2::DerivativeMapping* mapping) : impl(mapping) {fprintf(stderr, "Making Mapping\n");};
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
                      int nstokes, int nderiv, double* ssa,
                      double* total_extinction, double* emission_source,
                      double* f, double* leg_coeff, double* d_leg_coeff,
                      double* d_f, double* solar_irradiance);

    int get_derivative_mapping(const char* name,
                           DerivativeMapping** mapping);
};

struct Surface {
    std::unique_ptr<sasktran2::atmosphere::SurfaceInterface> impl;

    Surface(int nwavel, int nstokes);
};

struct Atmosphere {
    std::unique_ptr<sasktran2::atmosphere::AtmosphereInterface> impl;

    Atmosphere(AtmosphereStorage* storage, Surface* surface,
               bool calculate_derivatives);
};

struct OutputC {
    std::unique_ptr<sasktran2::OutputInterface> impl;

    OutputC(double* radiance, int nrad, int nstokes);
};

struct ViewingGeometry {
    sasktran2::viewinggeometry::ViewingGeometryContainer impl;

    ViewingGeometry() {}
};
