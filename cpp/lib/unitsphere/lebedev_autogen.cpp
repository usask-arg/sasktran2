#include <sasktran2/internal_common.h>
#include <sasktran2/math/unitsphere/lebedev.h>

#include <sasktran2/math/unitsphere/lebedev/sphere_6.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_14.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_26.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_38.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_50.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_74.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_86.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_110.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_146.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_170.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_194.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_230.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_266.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_302.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_350.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_434.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_590.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_770.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_974.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_1202.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_1454.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_1730.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_2030.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_2354.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_2702.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_3074.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_3470.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_3890.h>

namespace sasktran2::math::unitsphere::lebedev {
    void get_lebedev_data(int npoints, Eigen::MatrixXd& result,
                          Eigen::MatrixXi& faces) {
        if (npoints == 6) {
            result = Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_6, 4, 6);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_6, 3,
                                                      g_lebedev_num_faces_6);
        } else if (npoints == 14) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_14, 4, 14);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_14, 3,
                                                      g_lebedev_num_faces_14);
        } else if (npoints == 26) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_26, 4, 26);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_26, 3,
                                                      g_lebedev_num_faces_26);
        } else if (npoints == 38) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_38, 4, 38);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_38, 3,
                                                      g_lebedev_num_faces_38);
        } else if (npoints == 50) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_50, 4, 50);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_50, 3,
                                                      g_lebedev_num_faces_50);
        } else if (npoints == 74) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_74, 4, 74);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_74, 3,
                                                      g_lebedev_num_faces_74);
        } else if (npoints == 86) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_86, 4, 86);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_86, 3,
                                                      g_lebedev_num_faces_86);
        } else if (npoints == 110) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_110, 4, 110);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_110, 3,
                                                      g_lebedev_num_faces_110);
        } else if (npoints == 146) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_146, 4, 146);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_146, 3,
                                                      g_lebedev_num_faces_146);
        } else if (npoints == 170) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_170, 4, 170);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_170, 3,
                                                      g_lebedev_num_faces_170);
        } else if (npoints == 194) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_194, 4, 194);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_194, 3,
                                                      g_lebedev_num_faces_194);
        } else if (npoints == 230) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_230, 4, 230);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_230, 3,
                                                      g_lebedev_num_faces_230);
        } else if (npoints == 266) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_266, 4, 266);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_266, 3,
                                                      g_lebedev_num_faces_266);
        } else if (npoints == 302) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_302, 4, 302);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_302, 3,
                                                      g_lebedev_num_faces_302);
        } else if (npoints == 350) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_350, 4, 350);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_350, 3,
                                                      g_lebedev_num_faces_350);
        } else if (npoints == 434) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_434, 4, 434);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_434, 3,
                                                      g_lebedev_num_faces_434);
        } else if (npoints == 590) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_590, 4, 590);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_590, 3,
                                                      g_lebedev_num_faces_590);
        } else if (npoints == 770) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_770, 4, 770);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_770, 3,
                                                      g_lebedev_num_faces_770);
        } else if (npoints == 974) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_974, 4, 974);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_974, 3,
                                                      g_lebedev_num_faces_974);
        } else if (npoints == 1202) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_1202, 4, 1202);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_1202, 3,
                                                      g_lebedev_num_faces_1202);
        } else if (npoints == 1454) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_1454, 4, 1454);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_1454, 3,
                                                      g_lebedev_num_faces_1454);
        } else if (npoints == 1730) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_1730, 4, 1730);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_1730, 3,
                                                      g_lebedev_num_faces_1730);
        } else if (npoints == 2030) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_2030, 4, 2030);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_2030, 3,
                                                      g_lebedev_num_faces_2030);
        } else if (npoints == 2354) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_2354, 4, 2354);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_2354, 3,
                                                      g_lebedev_num_faces_2354);
        } else if (npoints == 2702) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_2702, 4, 2702);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_2702, 3,
                                                      g_lebedev_num_faces_2702);
        } else if (npoints == 3074) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_3074, 4, 3074);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_3074, 3,
                                                      g_lebedev_num_faces_3074);
        } else if (npoints == 3470) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_3470, 4, 3470);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_3470, 3,
                                                      g_lebedev_num_faces_3470);
        } else if (npoints == 3890) {
            result =
                Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_3890, 4, 3890);
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_3890, 3,
                                                      g_lebedev_num_faces_3890);
        } else {
            spdlog::error(
                "Requested number of Lebedev quadrature points does not exist");
            throw std::runtime_error(
                "Requested number of Lebedev quadrature points does not exist");
        }
    }
} // namespace sasktran2::math::unitsphere::lebedev
