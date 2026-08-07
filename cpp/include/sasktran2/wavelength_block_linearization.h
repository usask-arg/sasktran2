#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/wavelength_block.h>

namespace sasktran2 {
    /** Optional capability implemented by sources that linearize a complete
     * wavelength block together. Keeping this separate from the ordinary
     * source interface makes batched products an explicit, reviewable
     * extension rather than part of every source implementation. */
    class WavelengthBlockLinearizationInterface {
      public:
        virtual ~WavelengthBlockLinearizationInterface() = default;

        virtual bool select_calculated_block_lane(int, int) { return false; }

        virtual bool restore_forward_state_block(const WavelengthBlock<>&,
                                                 int) {
            return false;
        }

        virtual bool prepare_jvp_block(const WavelengthBlock<>&, int,
                                       const std::vector<Eigen::VectorXd>&) {
            return false;
        }

        /** Selects a zero-based lane from the block passed to
         * prepare_jvp_block. */
        virtual void select_jvp_block_lane(int, int) {}

        virtual bool begin_vjp_block(const WavelengthBlock<>&, int,
                                     Eigen::Ref<const Eigen::VectorXi>) {
            return false;
        }

        /** Selects a zero-based lane from the block passed to begin_vjp_block.
         */
        virtual void select_vjp_block_lane(int, int) {}
        /** Stages the cotangent for a zero-based lane in the active block. */
        virtual void stage_vjp_block_lane(int, int) {}

        virtual void finalize_vjp_block(const WavelengthBlock<>&, int,
                                        std::vector<Eigen::VectorXd>&) const {}
    };
} // namespace sasktran2
