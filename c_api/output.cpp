#include "output.h"
#include "sasktran2/config.h"
#include "internal_types.h"
#include <sasktran2.h>

OutputC::OutputC(double* radiance, int nrad, int nstokes) {
    Eigen::Map<Eigen::VectorXd> radiance_map(radiance, nrad);

    if (nstokes == 1) {
        impl = std::make_unique<sasktran2::OutputC<1>>(radiance_map);
    } else if (nstokes == 3) {
        impl = std::make_unique<sasktran2::OutputC<3>>(radiance_map);
    } else {
        // Handle error case
        impl = nullptr;
    }
}

extern "C" {
OutputC* sk_output_create(double* radiance, int nrad, int nstokes) {
    return new OutputC(radiance, nrad, nstokes);
}

void sk_output_destroy(OutputC* output) { delete output; }
}
