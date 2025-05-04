#include <limits>
#include <float.h>
#include <sasktran2/math/geodetic.h>
#include <sasktran2/math/trig.h>

static const double MACHINE_PRECISION = 100.0 * DBL_EPSILON;

namespace sasktran2::math::geodetic {
    class HeightOffsetEvaluator {
        Geodetic* m_geoid;              // The geoid used for calculations
        Eigen::Vector3d m_tangentpoint; // The location of the tangent point of
                                        // this look vector
        Eigen::Vector3d
            m_look; // The look direction as a unit vector (from the satellite).
        double m_H;

      public:
        HeightOffsetEvaluator(Geodetic* geoid, const Eigen::Vector3d& tanpoint,
                              const Eigen::Vector3d& look, double H);
        bool FindCrossingPoint(double lmin, double lmax,
                               Eigen::Vector3d* entrypoint);
        double operator()(double x);
    };

    Geodetic::Geodetic(double Re, double F)
        : m_geoid_Re(Re), m_geoid_F(F), m_is_valid(false),
          m_geodetic_altitude(std::numeric_limits<double>::quiet_NaN()),
          m_geodetic_latitude(std::numeric_limits<double>::quiet_NaN()),
          m_geodetic_longitude(std::numeric_limits<double>::quiet_NaN()) {
        m_local_south.setConstant(std::numeric_limits<double>::quiet_NaN());
        m_local_up.setConstant(std::numeric_limits<double>::quiet_NaN());
        m_local_west.setConstant(std::numeric_limits<double>::quiet_NaN());
        m_geocentric_location.setConstant(
            std::numeric_limits<double>::quiet_NaN());
    }

    Geodetic::~Geodetic() {}

    double Geodetic::altitude() const {
        if (!m_is_valid) {
            spdlog::warn("Accessing Geodetic::altitude when the geodetic has "
                         "not been initialized");
        }
        return m_geodetic_altitude;
    }

    double Geodetic::latitude() const {
        if (!m_is_valid) {
            spdlog::warn("Accessing Geodetic::latitude when the geodetic has "
                         "not been initialized");
        }
        return m_geodetic_latitude;
    }

    double Geodetic::longitude() const {
        if (!m_is_valid) {
            spdlog::warn("Accessing Geodetic::longitude when the geodetic has "
                         "not been initialized");
        }
        return m_geodetic_longitude;
    }

    std::pair<Eigen::Vector3d, Eigen::Vector3d>
    Geodetic::altitude_intercepts(double altitude,
                                  const Eigen::Vector3d& observer,
                                  const Eigen::Vector3d& look_vector) {
        Eigen::Vector3d entrypoint, exitpoint;

        Eigen::Vector3d tangentpoint; // The location of the tangent point of
                                      // this look vector
        double Ht;
        double Rmin;
        double Rmax;
        bool ok;
        bool ok1, ok2;
        double lmin, lmax;

        this->from_tangent_point(observer, look_vector);
        tangentpoint = location();
        Ht = this->altitude();
        ok = (altitude > Ht); // and make sure their is a valid solution
        if (ok)               // if there is then
        {
            if (m_geoid_F == 0) // If we have a pure sphere
            {                   // then the solution
                lmax = sqrt((2 * m_geoid_Re + altitude + Ht) *
                            (altitude - Ht)); // is considerably simplified
                Eigen::Vector3d d = look_vector.normalized();
                entrypoint = tangentpoint - lmax * d;
                exitpoint = tangentpoint + lmax * d;
            } else {
                HeightOffsetEvaluator evaluator(
                    this, tangentpoint, look_vector,
                    altitude); // Create the height offset evaluator object

                Rmin =
                    m_geoid_Re * (1.0 - m_geoid_F); // Get the Semi Minor axis
                Rmax = m_geoid_Re;                  // Get the Semi major axis
                lmin = 0.9 * sqrt((altitude - Ht) *
                                  (altitude + Ht +
                                   2 * Rmin)); // Minimum distance from tangent
                                               // point to shell height
                lmax = 1.1 * sqrt((altitude - Ht) *
                                  (altitude + Ht +
                                   2 * Rmax)); // Maximum distance from tangent
                                               // point to intersection

                ok1 = evaluator.FindCrossingPoint(
                    lmin, lmax, &entrypoint); // Find the shell from tangent
                                              // point towards observer
                ok2 = evaluator.FindCrossingPoint(
                    -lmin, -lmax, &exitpoint); // Find the shell from tangent
                                               // point away from observer
            }
            ok = ok1 && ok2;
            if (!ok) {
                spdlog::warn("Geodetic::altitude_intercepts, Error retrieving "
                             "entry and exit points for for height (%f)",
                             (double)altitude);
            }
        }
        if (!ok) {
            entrypoint.setZero();
            exitpoint.setZero();
        }
        return std::make_pair(entrypoint, exitpoint);
    }

    Eigen::Vector3d Geodetic::local_south() const {
        if (!m_is_valid) {
            spdlog::warn(
                "Accessing Geodetic::local_south when the geodetic has "
                "not been initialized");
        }
        return m_local_south;
    }

    Eigen::Vector3d Geodetic::local_up() const {
        if (!m_is_valid) {
            spdlog::warn("Accessing Geodetic::local_up when the geodetic has "
                         "not been initialized");
        }
        return m_local_up;
    }

    Eigen::Vector3d Geodetic::local_west() const {
        if (!m_is_valid) {
            spdlog::warn("Accessing Geodetic::local_west when the geodetic has "
                         "not been initialized");
        }
        return m_local_west;
    }

    Eigen::Vector3d Geodetic::location() const {
        if (!m_is_valid) {
            spdlog::warn("Accessing Geodetic::location when the geodetic has "
                         "not been initialized");
        }
        return m_geocentric_location;
    }

    void Geodetic::from_lat_lon_alt(double latitude, double longitude,
                                    double altitude) {
        m_geodetic_longitude = sasktran2::math::inrange(longitude, 360.0);
        m_geodetic_latitude = latitude;
        m_geodetic_altitude = altitude;

        double lat_rad =
            sasktran2::math::degrees_to_radians(m_geodetic_latitude);
        double lon_rad =
            sasktran2::math::degrees_to_radians(m_geodetic_longitude);

        double cosphi = cos(lat_rad);
        double sinphi = sin(lat_rad);
        double coslam = cos(lon_rad);
        double sinlam = sin(lon_rad);

        double fm1 = (1.0 - m_geoid_F);
        double fm12 = fm1 * fm1;

        double C = 1.0 / sqrt(cosphi * cosphi + fm12 * sinphi * sinphi);
        double S = fm12 * C;
        double P = (m_geoid_Re * C + m_geodetic_altitude) * cosphi;

        m_geocentric_location << P * coslam, P * sinlam,
            (m_geoid_Re * S + m_geodetic_altitude) * sinphi;

        update_local_coords();

        m_is_valid = true;
    }

    Eigen::Vector3d
    Geodetic::from_tangent_altitude(double altitude,
                                    const Eigen::Vector3d& observer,
                                    const Eigen::Vector3d& boresight) {
        Eigen::Vector3d
            xunit; // x unit vector in local vertical (upward) direction
        Eigen::Vector3d yunit; // y unit vector paralle to velocity vector
        double radius;
        double re_geocentric;
        Eigen::Vector3d offset;
        bool ok;
        int numtries;

        xunit = observer.normalized(); // Get unit vector to spacecraft

        yunit = (boresight - (boresight.dot(xunit) * xunit))
                    .normalized(); // Get unit vector perpendicular to r (to
                                   // define look plane)
        from_xyz(observer);
        get_osculating_spheroid(&re_geocentric, &offset);
        radius = (observer - offset).norm();

        double theta =
            (re_geocentric + altitude) /
            radius; // Get angle to tangent altitude on spherical Earth
        if (theta > 1.0) {
            spdlog::warn("nxGeodetic::FromTangentAltitude, Your requested "
                         "height is too high for this code. Tangent point is "
                         "above or equal to observers altitude");
            theta = 1.0;
        }
        theta = acos(theta);
        Eigen::Vector3d lookv =
            yunit * cos(theta) -
            xunit * sin(theta); // and initialize the look vector
        numtries = 0;
        do {
            from_tangent_point(observer, lookv); // calculate the tangent point
            double newh =
                m_geodetic_altitude; // get the geodetic height o fthis point
            double dh = (newh - altitude); // get the difference in height from
                                           // guess to geodetic
            ok = (fabs(dh) < 0.1);         // see if we have converged
            if (!ok)                       // we haven't
            {                              // so
                double dtheta =
                    0.8 * dh /
                    (radius * sin(theta)); // adjust the estimate of theta
                theta += dtheta;           // get a new theta
                lookv = yunit * cos(theta) -
                        xunit * sin(theta); // and unpdate the unit vectors
            }                               // and that is that
            numtries++;
        } while (!ok && numtries < 100); // repeat until converged;
        return lookv;
    }

    void Geodetic::from_tangent_point(const Eigen::Vector3d& observer,
                                      const Eigen::Vector3d& look_vector) {
        const int numiter = 5;

        Eigen::Vector3d loc;
        Eigen::Vector3d ell;
        // Stretch coordinates along the z-direction

        double f = m_geoid_F;

        ell = look_vector.normalized();

        for (int i = 0; i < numiter; ++i) {
            double transformFactor = 1.0 / (1.0 - f);
            Eigen::Vector3d obsTransf(observer.x(), observer.y(),
                                      observer.z() * transformFactor);
            Eigen::Vector3d ellTransf(look_vector.x(), look_vector.y(),
                                      look_vector.z() * transformFactor);

            // The radius is minimized at distance dot(look,obs) along the look
            // direction
            double s = -ellTransf.dot(obsTransf) / ellTransf.dot(ellTransf);

            // Set point, transforming back to ellipsoid
            loc = observer + s * ell;

            // Set as geodetic location, return to caller
            from_xyz(loc);

            f = m_geoid_F * (m_geoid_Re / (m_geoid_Re + m_geodetic_altitude));
        }
    }

    void Geodetic::get_osculating_spheroid(double* radius,
                                           Eigen::Vector3d* offset) {
        double oldheight = m_geodetic_altitude;
        double y0;
        double x0;
        double a;
        double b;
        double r;
        double theta;
        double dx;
        double dy;
        double a2y0;
        double b2x0;
        double a2;
        double b2;
        Eigen::Vector3d xunit(cosd(m_geodetic_longitude),
                              sind(m_geodetic_longitude), 0);
        Eigen::Vector3d yunit(0, 0, 1);

        from_lat_lon_alt(m_geodetic_latitude, m_geodetic_longitude, 0.0);

        y0 = m_geocentric_location.z(); // Get the vertical component of the
                                        // point on the surface of the geoid
        x0 = sqrt(sqr(m_geocentric_location.x()) +
                  sqr(m_geocentric_location
                          .y())); // Get the horizontal component of the
                                  // point on the surface of the geoid
        a = m_geoid_Re;
        b = m_geoid_Re * (1 - m_geoid_F);
        a2 = a * a;
        b2 = b * b;
        a2y0 = a2 * y0;
        b2x0 = b2 * x0;
        r = 1.0 / (a * b) *
            pow((a2y0 * y0 / b2 + b2x0 * x0 / a2),
                1.5); // Get the radius of curvature at this location
        theta = atan2(a2y0, b2x0); // Get the angle of the gradient of the
                                   // vertical at the surface of the geoid
        dx = r * cos(theta); // Get the X offset of the center of curvature from
                             // the surface point.
        dy = r * sin(theta); // Get the Y offset of the center of curvature from
                             // the surface point.
        *offset = m_geocentric_location - dy * yunit - dx * xunit;
        *radius = r;
        from_lat_lon_alt(m_geodetic_latitude, m_geodetic_longitude,
                         oldheight); // restore the original location
    }

    void Geodetic::from_xyz(const Eigen::Vector3d& location) {
        m_geocentric_location = location;

        double r;
        double x = m_geocentric_location.x();
        double y = m_geocentric_location.y();
        double z = m_geocentric_location.z();

        m_geodetic_longitude = sasktran2::math::inrange(atan2d(y, x), 360.0);
        if (m_geoid_F == 0) {
            r = sqrt(x * x + y * y + z * z);
            m_geodetic_latitude = sasktran2::math::asind(z / r);
            m_geodetic_altitude = r - m_geoid_Re;
        } else {
            r = sqrt(x * x + y * y);

            exact_geocentric_to_geodetic(r, z, &m_geodetic_latitude,
                                         &m_geodetic_altitude);
            m_geodetic_latitude =
                sasktran2::math::radians_to_degrees(m_geodetic_latitude);
        }

        update_local_coords();

        m_is_valid = true;
    }

    void Geodetic::update_local_coords() {
        Eigen::Vector3d horiz(m_geocentric_location.x(),
                              m_geocentric_location.y(), 0.0);
        horiz = horiz.normalized();

        Eigen::Vector3d vertical(0.0, 0.0, 1.0);

        m_local_up = horiz * sasktran2::math::cosd(m_geodetic_latitude) +
                     vertical * sasktran2::math::sind(m_geodetic_latitude);
        m_local_south = horiz * sasktran2::math::sind(m_geodetic_latitude) -
                        vertical * sasktran2::math::cosd(m_geodetic_latitude);
        m_local_west = m_local_south.cross(m_local_up);
    }

    void Geodetic::exact_geocentric_to_geodetic(double r, double z, double* fi,
                                                double* h) {
        double G[2];
        double DS[2];
        double T[4];
        int numsolutions;
        bool getallsolns = false;
        double bestt;
        double guesst = 0.0;
        double mindt;
        double dt;
        double zangle;
        double xangle;
        double l;
        bool gotsolution = false;
        bool zisnegative =
            (z < 0); // Keep track of whether z is negative or not.

        double s;
        double v;
        double a = m_geoid_Re; // approx 6378137.0 for Earth;
        double fr = m_geoid_F; // approx 1.0/298.257222101 for Earth
        double b = a - a * fr; // b is guaranteed same sign as z

        // ---- see if solution is close to the z axis as this needs special
        // attention

        z = fabs(z); // only deal inth angles in the first quadrant (zisnegative
                     // tracks the other stuff)
        l = sqrt(z * z + r * r);
        if (l < 1000.0 * std::numeric_limits<double>::epsilon()) {
            *fi = 0.0;
            *h = -m_geoid_Re;
        } else {
            zangle = (r / l);
            xangle = (z / l);
            if (zangle < 0.01)      // If we are within 1/100 of a radian
            {                       // then
                getallsolns = true; // get all of the solutions
                guesst = 0.0;       // t = 0 is the +ve z axis
                if (zangle <
                    MACHINE_PRECISION) // If we are really close to the Z axis
                {                      // within numerical precision
                    gotsolution =
                        true;       // then we can go straight to the solution
                    *fi = Pi / 2.0; // get the positive solution
                    *h = z - b;     // at the pole
                }
            }

            // ---- see if solution is close to the x axis as this needs special
            // attention

            if (xangle < 0.01)      // if we are within 1/100 of a radian
            {                       // then
                getallsolns = true; // get all of the solutions
                guesst = 1.0;       // t=1 is the +ve x axis.
                if (xangle <
                    MACHINE_PRECISION) // similarly if really close to the
                {                      // equator
                    gotsolution = true;
                    *fi = 0.0;  // then set the latitude
                    *h = r - a; // and get the height
                }
            }

            if (!gotsolution) { //  Find solution to: t**4 + 2*E*t**3 + 2*F*t -
                                //  1 = 0
                double sqrtd;
                double E = ((z + b) * b / a - a) / r;
                double F = ((z - b) * b / a + a) / r;
                double P = (E * F + 1.0) *
                           (4.0 / 3.0); // Calculate the intermediate terms
                double Q = (E * E - F * F) *
                           2.0; // as specified in the Borkowski paper.
                double D = P * P * P + Q * Q;

                if (D >= 0.0) // if (D > 0) as it usually is above ~47km from
                              // center of Earth
                {             // then we will normally want
                    sqrtd = sqrt(D);          // the cubic resolvent technique
                    s = std::cbrt(sqrtd + Q); // discussed in Borkowski paper to
                    v = P / s -
                        s; // to evaluate the difference of the two cube roots
                    v = -(Q + Q + v * v * v) /
                        (3.0 * P); // This is the formula given in the paper
                                   // (equation 20)
                    if (v * v > fabs(P)) // However if accuracy is not that good
                    { // which happens about 47km from center of Earth
                        v = std::cbrt(sqrtd - Q) -
                            std::cbrt(
                                sqrtd +
                                Q); // then use the direct method (equation 14a)
                    }               // which considerably improves the accuracy
                }                   // otherwise if D < 0
                else                // which happens close to the
                {                   // center of the Earth
                    double psqrt = sqrt(-P); // then
                    v = 2.0 * psqrt *
                        cos(acos(Q / (P * psqrt)) /
                            3.0); // use this alternate formula given by author
                                  // (equation 14b)
                }

                double esqrt = sqrt(E * E + v);
                double g = 0.5 * (E + esqrt); // Now get the default solution
                                              // for G (equation 13 in paper)
                double ds =
                    g * g +
                    (F - v * g) / (g + g - E); // Get the contents of square
                                               // root brackets in equation 12
                double ssqrt = sqrt(ds);       // and get the square root
                double t =
                    ssqrt - g; // solve equation 12 for the default value of t.

                if (getallsolns) // If we are close to the axes then we wantthe
                                 // solution
                { // that is also close to the axis, so we have to pick from
                  // solutions
                    G[0] = g;   // Get the 1st solution for G.
                    DS[0] = ds; // Get the square root term corresponding to +ve
                                // sqrt of G
                    G[1] = 0.5 * (E - esqrt); // Get the 2nd solution for G
                    DS[1] =
                        G[1] * G[1] +
                        (F - v * G[1]) /
                            (G[1] + G[1] - E); // Get the square root term
                                               // correpsondig to -ve sqrt of G

                    T[0] = t; // The first solution is the one we have already
                              // calculated
                    T[1] = -ssqrt - g; // the second is with negative of sqrt of
                                       // first value of ds

                    if (DS[1] >= 0) // See if we have the other two real
                                    // solutions or if they are complex.
                    {               // we do have them
                        T[2] = sqrt(DS[1]) - G[1];  // so get them
                        T[3] = -sqrt(DS[1]) - G[1]; // and
                        numsolutions =
                            4; // flag that we have all four real solutions
                    }          // otherwise
                    else       // we only have
                    {          // two real solutions
                        numsolutions = 2; // and that is that
                    }
                    bestt = T[0]; // Now hunt for the solution closes to our
                                  // "desired" guess estimate
                    mindt = fabs(T[0] -
                                 guesst); // initialize with the first default
                                          // solution (its often best)
                    for (int i = 1; i < numsolutions;
                         i++) // and then hunt through the list
                    {         // of other solutions
                        dt = fabs(T[i] - guesst); // to see if we have any that
                                                  // are even better
                        if (dt < mindt)           // if we do
                        {                         // then
                            bestt = T[i];         // copy it over
                            mindt = dt; // and set up the new minimum distance
                        }               // that is that
                    }                   // check all of the solutions
                    t = bestt;          // copy over the best solution into t
                }

                double divisor =
                    2.0 * b * t; // now use the best value of parameter t
                double numerator = (1.0 - t * t) * a; // to get
                if (fabs(divisor) != 0)
                    *fi = atan(numerator / divisor); // the geodetic latitude
                else
                    *fi = Pi / 2.0; // check for 90 degree stuff
                *h = (r - a * t) * cos(*fi) +
                     (z - b) * sin(*fi); // and get the geodetic altitude.
            }                            // and that is that
            if (zisnegative)
                *fi = -(*fi); // reverse the geodetic latitude if z was
                              // originally negative.
        }
    }

    template <class T>
    double zbrent(T func, double x1, double x2, double tol, int* status) {
        const int ITMAX = 500;
        const double EPS = 3.0e-8;
        int iter;
        double a = x1;
        double b = x2;
        double c = 0;
        double d = 0;
        double e = 0;
        double min1, min2;
        double fa = func(a);
        double fb = func(b);
        double fc, p, q, r, s, tol1, xm;

        *status = 0;
        if (fb * fa > 0.0)
            *status = 1;
        fc = fb;
        for (iter = 1; iter <= ITMAX; iter++) {
            if (fb * fc > 0.0) {
                c = a;
                fc = fa;
                e = d = b - a;
            }
            if (fabs(fc) < fabs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }
            tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
            xm = 0.5 * (c - b);
            if (fabs(xm) <= tol1 || fb == 0.0)
                return b;
            if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
                s = fb / fa;
                if (a == c) {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                } else {
                    q = fa / fc;
                    r = fb / fc;
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }
                if (p > 0.0)
                    q = -q;
                p = fabs(p);
                min1 = 3.0 * xm * q - fabs(tol1 * q);
                min2 = fabs(e * q);
                if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                    e = d;
                    d = p / q;
                } else {
                    d = xm;
                    e = d;
                }
            } else {
                d = xm;
                e = d;
            }
            a = b;
            fa = fb;
            if (fabs(d) > tol1)
                b += d;
            else
                b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
            fb = func(b);
        }
        *status = 1;
        return -9999.00; // nrerror("Maximum number of iterations exceeded in
                         // ZBRENT");
    }

    HeightOffsetEvaluator::HeightOffsetEvaluator(
        Geodetic* geoid, const Eigen::Vector3d& tanpoint,
        const Eigen::Vector3d& look, double H) {
        m_geoid = geoid;
        m_tangentpoint = tanpoint;
        m_look = look;
        m_H = H;
    }

    double HeightOffsetEvaluator::operator()(double x) {
        Eigen::Vector3d location;

        location = m_tangentpoint - m_look * x;
        m_geoid->from_xyz(location);
        return (m_geoid->altitude() - m_H);
    }

    bool HeightOffsetEvaluator::FindCrossingPoint(double lmin, double lmax,
                                                  Eigen::Vector3d* entrypoint) {
        double delta = 0.05 * lmin;
        double answer;
        int status;
        HeightOffsetEvaluator& evaluator = *this;

        answer =
            evaluator(lmin);  // Get the height offset closest to tangent point
        while (answer >= 0.0) // if
        {
            lmin -= delta;
            answer = evaluator(lmin);
        }

        answer = evaluator(lmax);
        while (answer <= 0.0) {
            lmax += 0.05 * lmax;
            answer = evaluator(lmax);
        }
        answer = zbrent(evaluator, lmin, lmax, 0.1,
                        &status); /* An estimate to the min location*/
        if (status == 0) {
            evaluator(answer);
            *entrypoint = m_geoid->location();
        } else {
            entrypoint->setZero();
        }
        return (status == 0);
    }

} // namespace sasktran2::math::geodetic
