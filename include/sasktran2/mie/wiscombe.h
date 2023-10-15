#include <sasktran2/mie/mie.h>


namespace sasktran2::mie {

    #ifdef INCLUDE_WISCOMBE_MIE

    class WiscombeMie : public MieBase {
        private:
            bool m_perfect;   /** Wiscombe's "perfect" i.e. infinite refractive index */
            double m_mincut;  /** Minimum imaginary refractive index befere truncating to 0 */
            int m_ipolzn;     /** wiscombe ipolzn variable */
            int m_nmom;       /** Number of moments to calculate */
        public:
            WiscombeMie();
            ~WiscombeMie() {};

            void calculate(const Eigen::VectorXd& size_param, const Eigen::VectorXcd& refractive_index, const Eigen::VectorXd& angles, MieOutput& output) override;
    };

    #endif
}


