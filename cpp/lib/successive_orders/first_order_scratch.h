#pragma once

#include <Eigen/Core>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

namespace sasktran2::successive_orders {

    /** Temporary compact-scalar products, shared only by one calling thread.
     *
     * Physical wavelength caches never live here. Backing vectors grow without
     * shrinking between engines; callers receive only their logical spans.
     * A lease must outlive all views and any source-thread parallel region.
     */
    class FirstOrderScratchLease {
      public:
        using View = Eigen::Map<Eigen::VectorXd>;

        FirstOrderScratchLease(int solar_threads, Eigen::Index solar_size,
                               Eigen::Index table_size)
            : m_slot(&thread_slot()), m_solar_threads(solar_threads),
              m_solar_size(solar_size), m_table_size(table_size) {
            if (solar_threads < 0 || solar_size < 0 || table_size < 0) {
                throw std::invalid_argument(
                    "Invalid compact first-order scratch dimensions");
            }
            if (m_slot->busy) {
                // A nested product must not grow or overwrite the outer lease.
                m_fallback = std::make_unique<Storage>();
                m_storage = m_fallback.get();
            } else {
                m_storage = &m_slot->storage;
                m_slot->busy = true;
                m_shared = true;
            }
            try {
                if (static_cast<int>(m_storage->solar.size()) < solar_threads) {
                    m_storage->solar.resize(solar_threads);
                }
                for (int thread = 0; thread < solar_threads; ++thread) {
                    grow(m_storage->solar[thread], solar_size);
                }
                grow(m_storage->table, table_size);
            } catch (...) {
                release();
                throw;
            }
        }

        ~FirstOrderScratchLease() { release(); }
        FirstOrderScratchLease(const FirstOrderScratchLease&) = delete;
        FirstOrderScratchLease&
        operator=(const FirstOrderScratchLease&) = delete;
        FirstOrderScratchLease(FirstOrderScratchLease&&) = delete;
        FirstOrderScratchLease& operator=(FirstOrderScratchLease&&) = delete;

        View solar(int thread) {
            if (thread < 0 || thread >= m_solar_threads) {
                throw std::out_of_range("Invalid first-order scratch thread");
            }
            return {m_storage->solar[thread].data(), m_solar_size};
        }
        View table() { return {m_storage->table.data(), m_table_size}; }

        // These distinct buffers may be grown after obtaining solar/table
        // views; their allocation cannot invalidate either of those views.
        void prepare_endpoints(Eigen::Index size) {
            if (size < 0) {
                throw std::invalid_argument(
                    "Invalid first-order scratch endpoint count");
            }
            grow(m_storage->endpoint_extinction, size);
            grow(m_storage->endpoint_albedo, size);
            m_endpoint_size = size;
        }
        View endpoint_extinction() {
            return {m_storage->endpoint_extinction.data(), m_endpoint_size};
        }
        View endpoint_albedo() {
            return {m_storage->endpoint_albedo.data(), m_endpoint_size};
        }

        bool shared() const { return m_shared; }
        std::uint64_t allocation_id() const { return m_storage->id; }
        std::size_t solar_bytes() const {
            std::size_t result = 0;
            for (const auto& values : m_storage->solar) {
                result += bytes(values);
            }
            return result;
        }
        std::size_t table_bytes() const { return bytes(m_storage->table); }
        std::size_t endpoint_bytes() const {
            return bytes(m_storage->endpoint_extinction) +
                   bytes(m_storage->endpoint_albedo);
        }
        std::size_t storage_bytes() const {
            return solar_bytes() + table_bytes() + endpoint_bytes();
        }

      private:
        static std::uint64_t next_id() {
            static std::atomic<std::uint64_t> id{0};
            return id.fetch_add(1, std::memory_order_relaxed) + 1;
        }
        struct Storage {
            const std::uint64_t id = next_id();
            std::vector<Eigen::VectorXd> solar;
            Eigen::VectorXd table;
            Eigen::VectorXd endpoint_extinction;
            Eigen::VectorXd endpoint_albedo;
        };
        struct ThreadSlot {
            Storage storage;
            bool busy = false;
        };

        static ThreadSlot& thread_slot() {
            // Retained until the calling OS thread exits, independently of
            // engine lifetime. Concurrent OS threads have distinct arenas.
            static thread_local ThreadSlot slot;
            return slot;
        }
        static void grow(Eigen::VectorXd& values, Eigen::Index size) {
            if (values.size() < size) {
                values.resize(size);
            }
        }
        static std::size_t bytes(const Eigen::VectorXd& values) {
            return static_cast<std::size_t>(values.size()) * sizeof(double);
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
        int m_solar_threads;
        Eigen::Index m_solar_size;
        Eigen::Index m_table_size;
        Eigen::Index m_endpoint_size = 0;
        bool m_shared = false;
    };

} // namespace sasktran2::successive_orders
