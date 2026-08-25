#pragma once

#include <Eigen/Core>
#include <Eigen/Cholesky>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace sasktran2::successive_orders {

    struct FixedPointSettings {
        int maximum_iterations = 50;
        double relative_tolerance = 1.0e-6;
        double absolute_tolerance = 1.0e-12;
        double damping = 1.0;
        int anderson_depth = 3;

        void validate() const {
            if (maximum_iterations < 0) {
                throw std::invalid_argument("successive-orders maximum "
                                            "iterations must be non-negative");
            }
            if (!std::isfinite(relative_tolerance) ||
                relative_tolerance < 0.0 ||
                !std::isfinite(absolute_tolerance) ||
                absolute_tolerance < 0.0) {
                throw std::invalid_argument("successive-orders tolerances must "
                                            "be finite and non-negative");
            }
            if (!std::isfinite(damping) || damping <= 0.0 || damping > 1.0) {
                throw std::invalid_argument(
                    "successive-orders damping must lie in (0, 1]");
            }
            if (anderson_depth < 0) {
                throw std::invalid_argument(
                    "successive-orders Anderson depth must be non-negative");
            }
        }

        bool convergence_enabled() const {
            return relative_tolerance > 0.0 || absolute_tolerance > 0.0;
        }
    };

    enum class FixedPointTermination { tolerance, maximum_iterations };

    struct FixedPointDiagnostics {
        static constexpr double material_warning_relative_tolerance = 1.0e-3;
        static constexpr double material_warning_absolute_tolerance = 1.0e-12;

        FixedPointTermination termination =
            FixedPointTermination::maximum_iterations;
        int iterations = 0;
        double residual_norm = std::numeric_limits<double>::infinity();
        double state_scale = 0.0;
        double convergence_threshold = 0.0;

        bool converged() const {
            return termination == FixedPointTermination::tolerance;
        }

        bool residual_exceeds(double relative_tolerance,
                              double absolute_tolerance) const {
            return residual_norm >
                   absolute_tolerance + relative_tolerance * state_scale;
        }

        double relative_residual() const {
            if (state_scale > 0.0) {
                return residual_norm / state_scale;
            }
            return residual_norm == 0.0
                       ? 0.0
                       : std::numeric_limits<double>::infinity();
        }

        bool warrants_convergence_warning() const {
            return !converged() &&
                   residual_exceeds(material_warning_relative_tolerance,
                                    material_warning_absolute_tolerance);
        }
    };

    /** Reusable Type-II Anderson workspace.
     *
     * Storage scales as O(n * depth), rather than retaining the complete
     * fixed-point history. The small normal equation is regularized relative
     * to its diagonal scale; if it is singular or produces a non-finite step,
     * the solver safely falls back to a damped Picard update.
     */
    class FixedPointWorkspace {
      public:
        void resize(Eigen::Index state_size, int depth) {
            if (state_size < 0 || depth < 0) {
                throw std::invalid_argument(
                    "invalid successive-orders solver workspace dimensions");
            }
            if (m_mapped.size() != state_size) {
                m_mapped.resize(state_size);
                m_residual.resize(state_size);
                m_previous_state.resize(state_size);
                m_previous_residual.resize(state_size);
            }
            if (m_delta_states.rows() != state_size ||
                m_delta_states.cols() != depth) {
                m_delta_states.resize(state_size, depth);
                m_delta_residuals.resize(state_size, depth);
            }
            if (m_gram.rows() != depth) {
                m_gram.resize(depth, depth);
                m_rhs.resize(depth);
                m_coefficients.resize(depth);
            }
            m_depth = depth;
            reset_history();
        }

        void reset_history() {
            m_history_size = 0;
            m_history_start = 0;
            m_have_previous = false;
        }

        Eigen::VectorXd& mapped() { return m_mapped; }
        Eigen::VectorXd& residual() { return m_residual; }

      private:
        friend class FixedPointSolver;

        Eigen::VectorXd m_mapped;
        Eigen::VectorXd m_residual;
        Eigen::VectorXd m_previous_state;
        Eigen::VectorXd m_previous_residual;
        Eigen::MatrixXd m_delta_states;
        Eigen::MatrixXd m_delta_residuals;
        Eigen::MatrixXd m_gram;
        Eigen::VectorXd m_rhs;
        Eigen::VectorXd m_coefficients;
        int m_depth = 0;
        int m_history_size = 0;
        int m_history_start = 0;
        bool m_have_previous = false;
    };

    class FixedPointSolver {
      public:
        /** Solves x = map(x).
         *
         * `map` must have the signature `void(const Eigen::VectorXd&,
         * Eigen::VectorXd&)`. It may reuse external scratch, but must support
         * distinct input/output vectors. When the iteration limit is reached,
         * the returned state is the last damped or Anderson iterate.
         */
        template <typename Map>
        static FixedPointDiagnostics solve(Eigen::VectorXd& state, Map&& map,
                                           const FixedPointSettings& settings,
                                           FixedPointWorkspace& workspace) {
            settings.validate();
            workspace.resize(state.size(), settings.anderson_depth);

            FixedPointDiagnostics diagnostics;
            if (settings.maximum_iterations == 0) {
                diagnostics.residual_norm = 0.0;
                return diagnostics;
            }

            for (int iteration = 0; iteration < settings.maximum_iterations;
                 ++iteration) {
                map(state, workspace.m_mapped);
                if (workspace.m_mapped.size() != state.size()) {
                    throw std::logic_error(
                        "successive-orders fixed-point map changed state size");
                }
                workspace.m_residual.noalias() = workspace.m_mapped - state;

                diagnostics.state_scale =
                    std::max(state.lpNorm<Eigen::Infinity>(),
                             workspace.m_mapped.lpNorm<Eigen::Infinity>());
                diagnostics.residual_norm =
                    workspace.m_residual.lpNorm<Eigen::Infinity>();
                diagnostics.convergence_threshold =
                    settings.absolute_tolerance +
                    settings.relative_tolerance * diagnostics.state_scale;
                diagnostics.iterations = iteration + 1;

                if (!std::isfinite(diagnostics.residual_norm) ||
                    !workspace.m_mapped.allFinite()) {
                    throw std::runtime_error(
                        "non-finite successive-orders fixed-point iterate");
                }

                if (settings.convergence_enabled() &&
                    diagnostics.residual_norm <=
                        diagnostics.convergence_threshold) {
                    state = workspace.m_mapped;
                    diagnostics.termination = FixedPointTermination::tolerance;
                    return diagnostics;
                }

                update_history(state, workspace.m_residual, workspace);
                if (!anderson_update(state, settings.damping, workspace)) {
                    state.noalias() += settings.damping * workspace.m_residual;
                }
            }

            return diagnostics;
        }

      private:
        static void update_history(const Eigen::VectorXd& state,
                                   const Eigen::VectorXd& residual,
                                   FixedPointWorkspace& workspace) {
            if (workspace.m_have_previous && workspace.m_depth > 0) {
                int slot;
                if (workspace.m_history_size < workspace.m_depth) {
                    slot =
                        (workspace.m_history_start + workspace.m_history_size) %
                        workspace.m_depth;
                    ++workspace.m_history_size;
                } else {
                    slot = workspace.m_history_start;
                    workspace.m_history_start =
                        (workspace.m_history_start + 1) % workspace.m_depth;
                }
                workspace.m_delta_states.col(slot).noalias() =
                    state - workspace.m_previous_state;
                workspace.m_delta_residuals.col(slot).noalias() =
                    residual - workspace.m_previous_residual;
            }
            workspace.m_previous_state = state;
            workspace.m_previous_residual = residual;
            workspace.m_have_previous = true;
        }

        static bool anderson_update(Eigen::VectorXd& state, double damping,
                                    FixedPointWorkspace& workspace) {
            const int history = workspace.m_history_size;
            if (history == 0) {
                return false;
            }

            for (int row = 0; row < history; ++row) {
                const int row_slot =
                    (workspace.m_history_start + row) % workspace.m_depth;
                const auto delta_residual_row =
                    workspace.m_delta_residuals.col(row_slot);
                workspace.m_rhs(row) =
                    delta_residual_row.dot(workspace.m_residual);
                for (int col = 0; col <= row; ++col) {
                    const int col_slot =
                        (workspace.m_history_start + col) % workspace.m_depth;
                    const double value = delta_residual_row.dot(
                        workspace.m_delta_residuals.col(col_slot));
                    workspace.m_gram(row, col) = value;
                    workspace.m_gram(col, row) = value;
                }
            }

            auto gram = workspace.m_gram.topLeftCorner(history, history);
            const double scale = std::max(1.0, gram.diagonal().maxCoeff());
            gram.diagonal().array() += 1.0e-14 * scale;
            workspace.m_coefficients.head(history) =
                gram.ldlt().solve(workspace.m_rhs.head(history));
            if (!workspace.m_coefficients.head(history).allFinite()) {
                return false;
            }

            workspace.m_mapped.noalias() =
                state + damping * workspace.m_residual;
            for (int column = 0; column < history; ++column) {
                const int slot =
                    (workspace.m_history_start + column) % workspace.m_depth;
                workspace.m_mapped.noalias() -=
                    workspace.m_coefficients(column) *
                    (workspace.m_delta_states.col(slot) +
                     damping * workspace.m_delta_residuals.col(slot));
            }
            if (!workspace.m_mapped.allFinite()) {
                return false;
            }
            state = workspace.m_mapped;
            return true;
        }
    };

} // namespace sasktran2::successive_orders
