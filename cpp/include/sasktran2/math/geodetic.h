#pragma once

#include "sasktran2/internal_common.h"

namespace sasktran2::math::geodetic {
    /**
     * @brief A class that represents a location within an ellipsoid coordinate
     * system such as WGS84
     *
     */
    class Geodetic {
      private:
        const double m_geoid_Re;
        const double m_geoid_F;

        bool m_is_valid;

        double m_geodetic_latitude;
        double m_geodetic_longitude;
        double m_geodetic_altitude;

        Eigen::Vector3d m_local_up;
        Eigen::Vector3d m_local_south;
        Eigen::Vector3d m_local_west;

        Eigen::Vector3d m_geocentric_location;

        /**
         * @brief Updates m_local_up, m_local_south, m_local_west after the
         * geodetic is set up
         *
         */
        void update_local_coords();

        /**
         * @brief Converts geocentric equatorial (r) and polar (z) components to
         *geodetic latitude (fi) and height (h) in [m]
         *
         * Code was taken from nxgeodetic.cxx on 2023-11-23, the docstring for
         *that function follows below
         *
         * Program to transform Cartesian to geodetic coordinates based
         *  on the exact solution (K.M. Borkowski,1989, Bulletin Geodesique. 63,
         *50-56) Input :  r, z = equatorial [m] and polar [m] components Output:
         *fi, h = geodetic coord's (latitude [rad], height [m])
         *
         *  This code came from http://www.astro.uni.torun.pl/~kb/geod.for
         *  Looks like it works well for the most part.
         *  I have tested the code in 1 km incremements from center of earth to
         *900 km altitude and it looks good.  I have also adjusted the handling
         *of Z and X axis as the original code picks "the wrong" solution in the
         *vicinity of Z and X axes (albeit a perfectly good solution).
         *
         *  There are some tolerances that have to built into the system and
         *these can be tricky to get just right.  The two tolerances are the
         *angle at which we consider being fairly close to the X or Z axis so we
         *choose the solution closest to the Z or X axis value. Curently this is
         *set to 0.01 radians or approx 1/2 a degree.
         *
         * @param r
         * @param z
         * @param fi
         * @param h
         */
        void exact_geocentric_to_geodetic(double r, double z, double* fi,
                                          double* h);

      public:
        /**
         * @brief Construct a new Geodetic object
         *
         * @param Re Equatorial radius of the geoid
         * @param F Flattening factor of the geoid
         */
        Geodetic(double Re, double F);
        ~Geodetic();

        /**
         * @brief The altitude of the point above the ellipsoid surface
         *
         * @return double
         */
        double altitude() const;

        /**
         * @brief The latitude of the point on the ellipsoid surface
         *
         * @return double
         */
        double latitude() const;

        /**
         * @brief The longitude of the point on the ellipsoid surface
         *
         * @return double
         */
        double longitude() const;

        /**
         * @brief A unit vector at the point that points in the local south
         * direction
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d local_south() const;

        /**
         * @brief A unit vector at the point that points perpindicular to the
         * ellipsoid
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d local_up() const;

        /**
         * @brief A unit vector at the point that points in the local west
         * direction
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d local_west() const;

        /**
         * @brief The point in (x, y, z) cartesian coordinates
         *
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d location() const;

        /**
         * For a given position (observer) and look vector (look_vector),
         * calculates the ellipsoidal intersectoin points that have an altitude
         * of altitude
         *
         * @param altitude Altitude of the desired intersection points [m]
         * @param observer Position of the observer
         * @param look_vector Local look vector
         * @return std::pair<Eigen::Vector3d, Eigen::Vector3d>
         */
        std::pair<Eigen::Vector3d, Eigen::Vector3d>
        altitude_intercepts(double altitude, const Eigen::Vector3d& observer,
                            const Eigen::Vector3d& look_vector);

        /**
         * @brief Initializes the geodetic from a give latitude, longitude, and
         * altitude
         *
         * @param latitude
         * @param longitude
         * @param altitude
         */
        void from_lat_lon_alt(double latitude, double longitude,
                              double altitude);

        /**
         * @brief Initializes the geodetic from a tangent altitude, observer
         * position, and boresight direction
         *
         * The boresight is scanned in the vertical direction until it results
         * in a tangent point with the requested altitude
         *
         * @param altitude Altitude of the desired tangent point [m]
         * @param observer Observer position
         * @param boresight Bore sight direction.  The resulting look vector
         * will be in the plane formed by the boresight and observer positions
         *
         * @return Eigen::Vector3d The resulting look vector
         */
        Eigen::Vector3d from_tangent_altitude(double altitude,
                                              const Eigen::Vector3d& observer,
                                              const Eigen::Vector3d& boresight);

        /**
         * @brief Constructs the geodetic from the tangent point formed by an
         * observer position and look vector
         *
         * @param observer The observer position
         * @param look_vector The look vector
         */
        void from_tangent_point(const Eigen::Vector3d& observer,
                                const Eigen::Vector3d& look_vector);

        /**
         * @brief Constructs the geodetic directly from an (x, y, z) position
         *
         * @param location
         */
        void from_xyz(const Eigen::Vector3d& location);

        /**
         * @brief True if the geodetic has been initialized
         *
         * @return true
         * @return false
         */
        bool is_valid() const { return m_is_valid; }

        /**
         * @brief Get the osculating spheroid parameters
         *
         * Code was taken from tangentpoint.cxx on 2012-11-23, the docstring for
         *that function follows below
         *
         *
         * Calculates the spheroid that best fits the ellipsoid surface (height
         *=0) at the current location. This algorithm best fits the spheroid in
         *a latitudinal direction (ie North South)
         *
         *  The radius of curvature is given by:-
         *
         *             [       (dy)^2 ] ^3/2
         *             [   1 + (--)   ]
         *             [       (dx)   ]
         *      R =  -----------------------
         *                d2y
         *                ---
         *                dx2
         *
         *  For an ellipse:
         *  R = 1/ab [ a^2y^2/b^2 + b^2x^2/a^2]^3/2
         *
         * @param radius
         * @param offset
         */
        void get_osculating_spheroid(double* radius, Eigen::Vector3d* offset);
    };
} // namespace sasktran2::math::geodetic
