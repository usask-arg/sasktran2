#include "output.h"
#include "sasktran2/config.h"
#include "internal_types.h"
#include <sasktran2.h>

OutputC::OutputC(double* radiance, int nrad) {
    Eigen::Map<Eigen::VectorXd> radiance_map(radiance, nrad);

    impl = std::make_unique<sasktran2::OutputC<1>>(radiance_map);
}

extern "C" {
OutputC* sk_output_create(double* radiance, int nrad) {
    return new OutputC(radiance, nrad);
}

void sk_output_destroy(OutputC* output) { delete output; }
}
