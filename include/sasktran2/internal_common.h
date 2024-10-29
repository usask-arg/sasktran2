#pragma once

#include <algorithm>
#include <numeric>
#include <mutex>
#include <queue>
#include <math.h>
#include <complex>
#include <iostream>
#include <optional>

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
#define lapack_int int
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

#include <Eigen/src/misc/lapacke.h>

#endif
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/CXX11/Tensor>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>
