#pragma once

#include <algorithm>
#include <numeric>
#include <mutex>
#include <queue>
#include <math.h>
#include <complex>
#include <iostream>
#include <optional>

#ifdef SKTRAN_TRACY
#include <tracy/Tracy.hpp>
#else
#include <sasktran2/tracy_dummy.h>
#endif

#include <sasktran2/types.h>

// Minimimum size shells have to be to have contributions
#define MINIMUM_SHELL_SIZE_M 0.0001

// If dot(viewing, normal) < NADIR_VIEWING_CUTOFF then the ray is considered to
// be strictly nadir viewing and the refraction is not applied
#define NADIR_VIEWING_CUTOFF 0.999999

// Setup dependening on what linear algebra package is being linked
#ifdef SKTRAN_USE_MKL
// Using MKL for linear algebdra
#define EIGEN_USE_MKL_ALL 1

#include <mkl_lapacke.h>
#else
// Unsure if this is faster or not
// #define EIGEN_USE_BLAS 1
#ifdef SKTRAN_USE_ACCELERATE
// Using apple Accelerate for linear algebra, which doesn't have a LAPACKE
// interface
#ifdef ACCELERATE_LAPACK_ILP64
#define lapack_int long
#else
#define lapack_int int
#endif
// #include <vecLib.h>
#include <clapack.h>
#include <cblas.h>
#else
#define LAPACK_DISABLE_NAN_CHECK
// Using a standard LAPACKE compatible package
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

// If you set this then various asserts will be checked during runtime and
// error/warnings emitted, useful for debugging in Python or for calculations
// that are too slow in debug mode #define SASKTRAN_DEBUG_ASSERTS

#ifdef SKTRAN_SCIPY_BLAS

#define __EMSCRIPTEN__
#define dgbsv_ scipy_dgbsv_64_
#define dgbtrs_ scipy_dgbtrs_64_
#define dgeev_ scipy_dgeev_64_
#include <lapack.h>
#undef __EMSCRIPTEN__
#else

#include <Eigen/src/misc/lapacke.h>

#endif
#endif
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/CXX11/Tensor>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>
