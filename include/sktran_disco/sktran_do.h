#pragma once

// STL dependencies

#include <sasktran2/internal_common.h>

// Eigen settings

// Whether or not to use the Eigensolver from Eigen instead of dgeev call
// directly I think the Eigen EigenSolver is faster for small matrices but
// slower for large ones, and sometimes the Eigen EigenSolver can have precision
// problems so right now we just disable it
#define SASKTRAN_DISCO_USE_EIGEN_EIGENSOLVER true

// If set, then the pentadiagonal solver is used for NSTR=2, NSTOKES=1, only
// enabled if SKTRANDO_FULL_COMPILE is set
#define SASKTRAN_DISCO_ENABLE_PENTADIAGONAL true

// Debug settings to enable or di
#define SASKTRAN_DISCO_ENABLE_FULL_BACKPROP true

// SKTRAN_DO is templated over two main parameters, the first is NSTOKES which
// we always explicitly instantiate every class over.  DO is also templated over
// the number of streams, the default is -1 which is "dynamic" and allows for
// any number of streams.  We also can specialize specific values of the number
// of streams but this takes a long time to compile so we typically only do it
// on release.  Useful values for speed are 2, 4, and 16

// #define SASKTRAN_DISCO_FULL_COMPILE

#ifdef SASKTRAN_DISCO_FULL_COMPILE
#define SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(classname)                         \
    template class classname<1>;                                               \
    template class classname<3>;                                               \
    template class classname<1, 2>;                                            \
    template class classname<3, 2>;                                            \
    template class classname<1, 4>;                                            \
    template class classname<3, 4>;                                            \
    template class classname<1, 16>;                                           \
    template class classname<3, 16>;

#define SASKTRAN_DISCO_INSTANTIATE_TEMPLATE_STRUCT(classname)                  \
    template struct classname<1>;                                              \
    template struct classname<3>;                                              \
    template struct classname<1, 2>;                                           \
    template struct classname<3, 2>;                                           \
    template struct classname<1, 4>;                                           \
    template struct classname<3, 4>;                                           \
    template struct classname<1, 16>;                                          \
    template struct classname<3, 16>;

#else
#define SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(classname)                         \
    template class classname<1>;                                               \
    template class classname<3>;

#define SASKTRAN_DISCO_INSTANTIATE_TEMPLATE_STRUCT(classname)                  \
    template struct classname<1>;                                              \
    template struct classname<3>;
#endif

// SASKTRAN dependencies

// SASKTRAN-DO internal dependencies
#include "sktran_disco/sktran_do_linearization_types.h"
#include "sktran_disco/sktran_do_polarization_types.h"
#include "sktran_disco/sktran_do_types.h"
#include "sktran_disco/sktran_do_linalg.h"
#include "sktran_disco/sktran_do_memory.h"
#include "sktran_disco/sktran_do_quadrature.h"
#include "sktran_disco/sktran_do_specs.h"
#include "sktran_disco/sktran_do_surface.h"
#include "sktran_disco/sktran_do_testing.h"
#include "sktran_disco/sktran_do_properties.h"
#include "sktran_disco/sktran_do_pconfig.h"
#include "sktran_disco/sktran_do_opticallayer.h"
#include "sktran_disco/sktran_do_geometrylayerarray.h"
#include "sktran_disco/sktran_do_postprocessing.h"
#include "sktran_disco/sktran_do_layerarray.h"
#include "sktran_disco/sktran_do_rte.h"
