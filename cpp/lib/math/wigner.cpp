#include "sasktran2/math/wigner.h"
#include "sasktran2/rust.h"

namespace sasktran2::math {
       std::unique_ptr<WignerDCalculator> new_wigner_d_calculator(int m, int n) {

        std::cout << test() << std::endl;

        return std::make_unique<WignerDCalculator>(m, n);
    } 
}