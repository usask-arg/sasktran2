#pragma once

#include "problem.h"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>

namespace sasktran2::successive_orders {

    /** Ephemeral scalar solver scratch owned by one calling OS thread.
     *
     * These are the original owning Eigen vectors/matrices, with their exact
     * original logical shapes and strides. Problem/scattering prepare methods
     * resize them when required, and each fixed-point solve resets history.
     * Physical operators, wavelength states/forcing and public ray products
     * never live here. All references must finish before this lease ends.
     */
    class ScalarProblemWorkspaceLease {
      public:
        ScalarProblemWorkspaceLease() : m_slot(&thread_slot()) {
            if (m_slot->busy) {
                // Reentrant calls cannot resize or overwrite the outer solve.
                m_fallback = std::make_unique<Storage>();
                m_storage = m_fallback.get();
            } else {
                m_storage = &m_slot->storage;
                m_slot->busy = true;
                m_shared = true;
            }
        }
        ~ScalarProblemWorkspaceLease() { release(); }
        ScalarProblemWorkspaceLease(const ScalarProblemWorkspaceLease&) =
            delete;
        ScalarProblemWorkspaceLease&
        operator=(const ScalarProblemWorkspaceLease&) = delete;
        ScalarProblemWorkspaceLease(ScalarProblemWorkspaceLease&&) = delete;
        ScalarProblemWorkspaceLease&
        operator=(ScalarProblemWorkspaceLease&&) = delete;

        ProblemWorkspace<1>& workspace() { return m_storage->workspace; }
        const ProblemWorkspace<1>& workspace() const {
            return m_storage->workspace;
        }
        bool shared() const { return m_shared; }
        std::uint64_t allocation_id() const { return m_storage->id; }
        std::size_t problem_bytes() const {
            return workspace().storage_bytes();
        }
        std::size_t fixed_point_bytes() const {
            return workspace().fixed_point.storage_bytes();
        }
        std::size_t storage_bytes() const {
            return problem_bytes() + fixed_point_bytes();
        }

      private:
        static std::uint64_t next_id() {
            static std::atomic<std::uint64_t> id{0};
            return id.fetch_add(1, std::memory_order_relaxed) + 1;
        }
        struct Storage {
            const std::uint64_t id = next_id();
            ProblemWorkspace<1> workspace;
        };
        struct ThreadSlot {
            Storage storage;
            bool busy = false;
        };
        static ThreadSlot& thread_slot() {
            // Scratch survives individual engines but is freed at OS-thread
            // exit. Concurrent calling threads own independent workspaces.
            static thread_local ThreadSlot slot;
            return slot;
        }
        void release() noexcept {
            if (m_shared) {
                m_slot->busy = false;
                m_shared = false;
            }
        }

        ThreadSlot* m_slot;
        std::unique_ptr<Storage> m_fallback;
        Storage* m_storage = nullptr;
        bool m_shared = false;
    };

} // namespace sasktran2::successive_orders
