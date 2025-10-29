#pragma once

#include <sasktran2/internal_common.h>

namespace sasktran2::math::unitsphere::lebedev {
    // Packed as x0, y0, z0, w0, x1, y1, z1, w1, ...
    const static double g_lebedev_xyzw_6[4 * 6] = {
        1.0,  0.0,  0.0,  0.16666666666666671,
        -1.0, 0.0,  0.0,  0.16666666666666671,
        0.0,  1.0,  0.0,  0.16666666666666671,
        0.0,  -1.0, 0.0,  0.16666666666666671,
        0.0,  0.0,  1.0,  0.16666666666666671,
        0.0,  0.0,  -1.0, 0.16666666666666671};

    // Triangulated faces (CCW, outward), 3 * g_lebedev_num_faces_6 entries
    const static int g_lebedev_num_faces_6 = 8;
    const static int g_lebedev_faces_6[3 * g_lebedev_num_faces_6] = {
        5, 3, 1, 5, 0, 3, 4, 1, 3, 4, 3, 0, 2, 0, 5, 2, 5, 1, 2, 4, 0, 2, 1, 4};
} // namespace sasktran2::math::unitsphere::lebedev
