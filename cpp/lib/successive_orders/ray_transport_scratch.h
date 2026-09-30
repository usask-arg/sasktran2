#pragma once

#include "ray_transport.h"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>

namespace sasktran2::successive_orders {

    /** Scratch for a synchronous scalar line-of-sight derivative operation.
     *
     * JVP and VJP use the same owning value buffer at disjoint times. The VJP
     * layer workspace is also transient. Neither physical transport values nor
     * the resulting state/ray products live here. References must finish before
     * the lease ends. Concurrent OS threads and reentrant calls are isolated.
     */
    class ScalarRayTransportWorkspaceLease {
      public:
        ScalarRayTransportWorkspaceLease() : m_slot(&thread_slot()) {
            if (m_slot->busy) {
                m_fallback = std::make_unique<Storage>();
                m_storage = m_fallback.get();
            } else {
                m_storage = &m_slot->storage;
                m_slot->busy = true;
                m_shared = true;
            }
        }
        ~ScalarRayTransportWorkspaceLease() {
            if (m_shared) {
                m_slot->busy = false;
            }
        }
        ScalarRayTransportWorkspaceLease(
            const ScalarRayTransportWorkspaceLease&) = delete;
        ScalarRayTransportWorkspaceLease&
        operator=(const ScalarRayTransportWorkspaceLease&) = delete;
        ScalarRayTransportWorkspaceLease(ScalarRayTransportWorkspaceLease&&) =
            delete;
        ScalarRayTransportWorkspaceLease&
        operator=(ScalarRayTransportWorkspaceLease&&) = delete;

        Eigen::VectorXd& values() { return m_storage->values; }
        RayTransportWorkspace& workspace() { return m_storage->workspace; }
        bool shared() const { return m_shared; }
        std::uint64_t allocation_id() const { return m_storage->id; }
        std::size_t value_bytes() const {
            return static_cast<std::size_t>(m_storage->values.size()) *
                   sizeof(double);
        }
        std::size_t workspace_bytes() const {
            return m_storage->workspace.storage_bytes();
        }
        std::size_t storage_bytes() const {
            return value_bytes() + workspace_bytes();
        }

      private:
        static std::uint64_t next_id() {
            static std::atomic<std::uint64_t> id{0};
            return id.fetch_add(1, std::memory_order_relaxed) + 1;
        }
        struct Storage {
            const std::uint64_t id = next_id();
            Eigen::VectorXd values;
            RayTransportWorkspace workspace;
        };
        struct ThreadSlot {
            Storage storage;
            bool busy = false;
        };
        static ThreadSlot& thread_slot() {
            // Freed at OS-thread exit, independently of individual engines.
            static thread_local ThreadSlot slot;
            return slot;
        }

        ThreadSlot* m_slot;
        std::unique_ptr<Storage> m_fallback;
        Storage* m_storage = nullptr;
        bool m_shared = false;
    };

} // namespace sasktran2::successive_orders
