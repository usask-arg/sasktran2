#include <sasktran2/mie/mie.h>
#include <sasktran2/mie/wiscombe.h>

#if !defined(__WIN32__)
 	#define MIEV0 miev0_
#endif

extern "C" void  MIEV0  (	double*					XX,				// REAL
							std::complex<double>*	CREFIN,			// COMPLEX
							char*					PERFCT,			// LOGICAL
							double*					MIMCUT,			// REAL
							char*					ANYANG,			// LOGICAL
							int*					NUMANG,			// INTEGER
							double*					XMU,			// REAL(*)
							int*					NMOM,			// INTEGER
							int*					IPOLZN,			// INTEGER
							int*					MOMDIM,			// INTEGER
							char*					PRNT,			// LOGICAL(*)
							double*					QEXT,			// REAL
							double*					QSCA,			// REAL
							double*					GQSC,			// REAL
							double*					PMOM,			// REAL(0:MOMDIM,*)
							std::complex<double>*		SFORW,			// COMPLEX
							std::complex<double>*		SBACK,			// COMPLEX
                            std::complex<double>*		S1,				// COMPLEX (NUMANG))
                            std::complex<double>*		S2,				// COMPLEX (NUMANG)
                            std::complex<double>*		TFORW,			// COMPLEX (2)
                            std::complex<double>*		TBACK,			// COMPLEX (2)
							double*					SPIKE );		// REAL

namespace sasktran2::mie {
    WiscombeMie::WiscombeMie() :
        m_perfect(false),
        m_mincut(0),
        m_ipolzn(4321),
        m_nmom(0)
    {

    }
    void WiscombeMie::calculate(const Eigen::VectorXd& size_param, const Eigen::VectorXcd& refractive_index, const Eigen::VectorXd& angles, MieOutput& output) {
        output.values.resize(size_param.size(), angles.size());

        output.size_param = size_param;
        output.refractive_index = refractive_index;
        output.angles = angles;

        char	PRNT[8]    = {0,0,0,0,0,0,0,0};
        char	PERFECT[4];
        char	ANYANG[4];
        int		NUMANG     = (int)angles.size();
        int		MOMDIM	   = m_nmom + 1;
        char	boolval;

        boolval = m_perfect ? -1 : 0;
        memset( PERFECT, boolval, sizeof(PERFECT));
        boolval = true ? -1 : 0;
        memset( ANYANG, boolval, sizeof(ANYANG));

        double spike;

        Eigen::MatrixXd phase_moments(m_nmom+2, 4);

        for(int i = 0; i < size_param.size(); ++i) {
            double xx = size_param(i);
            std::complex<double> crefin = refractive_index(i);

            MIEV0(
                &xx,
                &crefin,
                PERFECT,
                &m_mincut,
                ANYANG,
                &NUMANG,
                const_cast<double*>(angles.data()),
                &m_nmom,
                &m_ipolzn,
                &MOMDIM,
                PRNT,
                &output.values.Qext(i),
                &output.values.Qsca(i),
                &output.values.GQsc(i),
                phase_moments.data(),
                &output.values.SFORW(i),
                &output.values.SBACK(i),
                output.values.S1.row(i).data(),
                output.values.S2.row(i).data(),
                output.values.TFORW.row(i).data(),
                output.values.TBACK.row(i).data(),
                &spike
            );

        }


    }
}