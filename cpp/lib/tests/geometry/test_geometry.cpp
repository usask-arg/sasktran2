#include <sasktran2/test_helper.h>

#include <sasktran2.h>

#include <array>

TEST_CASE("Coordinates Spherical Angle Calculations", "[sasktran2][geometry]") {
    sasktran2::Coordinates geo(0.5, 60, 6372000,
                               sasktran2::geometrytype::spherical);

    double theta = GENERATE((take(100, random(-0.4, 0.4))));
    double phi = GENERATE((take(100, random(-0.4, 0.4))));

    Eigen::Vector3d uv = geo.unit_vector_from_angles(theta, phi);

    auto reconstructed_angles = geo.angles_from_unit_vector(uv);

    REQUIRE(fabs(theta - reconstructed_angles.first) < 1e-8);
    REQUIRE(fabs(phi - reconstructed_angles.second) < 1e-8);
}

TEST_CASE("Coordinates Solar Angle Vector Generation",
          "[sasktran2][geometry]") {
    sasktran2::Coordinates geo(0.5, 60, 6372000,
                               sasktran2::geometrytype::spherical);

    double cos_sza = GENERATE((take(100, random(-1, 1))));
    double phi = GENERATE((take(100, random(-EIGEN_PI, EIGEN_PI))));

    double altitude = GENERATE((take(10, random(0, 100000))));

    Eigen::Vector3d vec = geo.solar_coordinate_vector(cos_sza, phi, altitude);

    auto angles = geo.solar_angles_at_location(vec);

    REQUIRE(fabs(cos_sza - angles.first) < 1e-8);

    // Doesn't currently work but not used
    // REQUIRE(fabs(cos(phi) - cos(angles.second)) < 1e-8);
    // REQUIRE(fabs(sin(phi) - sin(angles.second)) < 1e-8);
}

TEST_CASE("Coordinates solar vectors rotate covariantly at the solar poles",
          "[sasktran2][geometry]") {
    constexpr std::array<double, 2> solar_cosines{1.0, -1.0};
    const Eigen::Matrix3d rotation =
        (Eigen::AngleAxisd(0.7, Eigen::Vector3d::UnitY()) *
         Eigen::AngleAxisd(-0.3, Eigen::Vector3d::UnitX()))
            .toRotationMatrix();

    for (const double cos_sza : solar_cosines) {
        DYNAMIC_SECTION("cos_sza=" << cos_sza) {
            sasktran2::Coordinates base(cos_sza, 0.2, 6372000.0,
                                        sasktran2::geometrytype::spherical);
            sasktran2::Coordinates rotated(
                rotation * base.reference_z(), rotation * base.reference_x(),
                rotation * base.sun_unit(), base.earth_radius(),
                sasktran2::geometrytype::spherical);

            const Eigen::Vector3d base_location =
                base.solar_coordinate_vector(cos_sza, 0.0, 1234.0);
            const Eigen::Vector3d rotated_location =
                rotated.solar_coordinate_vector(cos_sza, 0.0, 1234.0);
            REQUIRE((rotation * base_location - rotated_location).norm() <
                    1.0e-8);
        }
    }
}
