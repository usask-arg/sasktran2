#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace sasktran2::successive_orders {

    /** A finalized CSR slot and the unchanged representation of its weight. */
    template <typename Index> class CompactSourceWeight {
      public:
        CompactSourceWeight(std::uint32_t index, double weight) {
            static_assert(std::is_unsigned_v<Index>);
            if (index > std::numeric_limits<Index>::max()) {
                throw std::out_of_range("Source interpolation slot overflow");
            }
            m_index = static_cast<Index>(index);
            std::memcpy(m_weight.data(), &weight, sizeof(weight));
        }
        std::uint32_t row_inner_index() const { return m_index; }
        double weight() const {
            double value;
            std::memcpy(&value, m_weight.data(), sizeof(value));
            return value;
        }

      private:
        Index m_index;
        std::array<std::byte, sizeof(double)> m_weight;
    };

    static_assert(sizeof(CompactSourceWeight<std::uint8_t>) == 9);
    static_assert(sizeof(CompactSourceWeight<std::uint16_t>) == 10);

    template <typename T> class SourceWeightArrayView {
      public:
        SourceWeightArrayView(const T* data, std::size_t size)
            : m_data(size == 0 ? nullptr : data), m_size(size) {}
        const T* data() const { return m_data; }
        const T* begin() const { return m_data; }
        const T* end() const { return m_size == 0 ? m_data : m_data + m_size; }
        std::size_t size() const { return m_size; }
        bool empty() const { return m_size == 0; }
        const T& operator[](std::size_t index) const {
            if (index >= m_size) {
                throw std::out_of_range("Source interpolation view index");
            }
            return m_data[index];
        }

      private:
        const T* m_data;
        std::size_t m_size;
    };

    /** Immutable view of any supported slot width.
     *
     * Hot kernels visit a typed range once, so width selection does not occur
     * for each multiply/add. Legacy iteration returns values rather than
     * references into a storage format with a different record size.
     */
    template <typename WideWeight> class SourceWeightView {
      public:
        using ByteWeight = CompactSourceWeight<std::uint8_t>;
        using ShortWeight = CompactSourceWeight<std::uint16_t>;

        class Iterator {
          public:
            using iterator_category = std::input_iterator_tag;
            using value_type = WideWeight;
            using difference_type = std::ptrdiff_t;
            using reference = WideWeight;
            using pointer = void;

            WideWeight operator*() const {
                if (m_index >= m_size) {
                    throw std::out_of_range("Source interpolation iterator");
                }
                if (m_index_bytes == 1) {
                    const auto& value =
                        static_cast<const ByteWeight*>(m_data)[m_index];
                    return {static_cast<int>(value.row_inner_index()),
                            value.weight()};
                }
                if (m_index_bytes == 2) {
                    const auto& value =
                        static_cast<const ShortWeight*>(m_data)[m_index];
                    return {static_cast<int>(value.row_inner_index()),
                            value.weight()};
                }
                return static_cast<const WideWeight*>(m_data)[m_index];
            }
            Iterator& operator++() {
                ++m_index;
                return *this;
            }
            Iterator operator++(int) {
                auto previous = *this;
                ++*this;
                return previous;
            }
            bool operator==(const Iterator& other) const {
                return m_data == other.m_data && m_index == other.m_index &&
                       m_index_bytes == other.m_index_bytes;
            }
            bool operator!=(const Iterator& other) const {
                return !(*this == other);
            }

          private:
            friend class SourceWeightView;
            Iterator(const void* data, std::size_t size, std::size_t index,
                     unsigned char index_bytes)
                : m_data(data), m_size(size), m_index(index),
                  m_index_bytes(index_bytes) {}
            const void* m_data;
            std::size_t m_size;
            std::size_t m_index;
            unsigned char m_index_bytes;
        };

        SourceWeightView() = default;
        template <typename T>
        SourceWeightView(const std::vector<T>& values, std::size_t offset,
                         std::size_t size)
            : m_size(size) {
            static_assert(std::is_same_v<T, WideWeight> ||
                          std::is_same_v<T, ByteWeight> ||
                          std::is_same_v<T, ShortWeight>);
            if (offset > values.size() || size > values.size() - offset) {
                throw std::out_of_range("Source interpolation view bounds");
            }
            m_data = size == 0 ? nullptr : values.data() + offset;
            if constexpr (std::is_same_v<T, ByteWeight>) {
                m_index_bytes = 1;
            } else if constexpr (std::is_same_v<T, ShortWeight>) {
                m_index_bytes = 2;
            }
        }
        std::size_t size() const { return m_size; }
        bool empty() const { return m_size == 0; }
        unsigned char index_bytes() const { return m_index_bytes; }
        Iterator begin() const { return {m_data, m_size, 0, m_index_bytes}; }
        Iterator end() const { return {m_data, m_size, m_size, m_index_bytes}; }
        WideWeight operator[](std::size_t index) const {
            return *Iterator(m_data, m_size, index, m_index_bytes);
        }
        template <typename Function>
        decltype(auto) visit(Function&& function) const {
            if (m_index_bytes == 1) {
                return std::forward<Function>(function)(
                    SourceWeightArrayView<ByteWeight>(
                        static_cast<const ByteWeight*>(m_data), m_size));
            }
            if (m_index_bytes == 2) {
                return std::forward<Function>(function)(
                    SourceWeightArrayView<ShortWeight>(
                        static_cast<const ShortWeight*>(m_data), m_size));
            }
            return std::forward<Function>(function)(
                SourceWeightArrayView<WideWeight>(
                    static_cast<const WideWeight*>(m_data), m_size));
        }

      private:
        const void* m_data = nullptr;
        std::size_t m_size = 0;
        unsigned char m_index_bytes = 4;
    };

    /** Wide during construction, then the smallest exact finalized CSR slot.
     * Rows too large for the compact formats retain the wide representation.
     */
    template <typename WideWeight> class SourceWeightStorage {
      public:
        using value_type = WideWeight;
        using View = SourceWeightView<WideWeight>;
        using ByteWeight = typename View::ByteWeight;
        using ShortWeight = typename View::ShortWeight;

        SourceWeightStorage() = default;
        SourceWeightStorage(std::initializer_list<WideWeight> values)
            : m_values(std::vector<WideWeight>(values)) {}
        SourceWeightStorage&
        operator=(std::initializer_list<WideWeight> values) {
            m_values = std::vector<WideWeight>(values);
            return *this;
        }
        std::vector<WideWeight>& wide_values() {
            auto* values = std::get_if<std::vector<WideWeight>>(&m_values);
            if (values == nullptr) {
                throw std::logic_error("Finalized source storage is immutable");
            }
            return *values;
        }
        const std::vector<WideWeight>& wide_values() const {
            const auto* values =
                std::get_if<std::vector<WideWeight>>(&m_values);
            if (values == nullptr) {
                throw std::logic_error(
                    "Source storage no longer has wide records");
            }
            return *values;
        }
        void reserve(std::size_t size) { wide_values().reserve(size); }
        void append(const std::vector<WideWeight>& values) {
            auto& wide = wide_values();
            wide.insert(wide.end(), values.begin(), values.end());
        }
        std::size_t size() const {
            return std::visit([](const auto& values) { return values.size(); },
                              m_values);
        }
        bool empty() const { return size() == 0; }
        std::size_t capacity() const {
            return std::visit(
                [](const auto& values) { return values.capacity(); }, m_values);
        }
        std::size_t capacity_bytes() const {
            return std::visit(
                [](const auto& values) {
                    using T =
                        typename std::decay_t<decltype(values)>::value_type;
                    return values.capacity() * sizeof(T);
                },
                m_values);
        }
        std::size_t element_bytes() const {
            return std::visit(
                [](const auto& values) {
                    using T =
                        typename std::decay_t<decltype(values)>::value_type;
                    return sizeof(T);
                },
                m_values);
        }
        bool is_compact() const { return m_values.index() != 0; }
        void shrink_to_fit() {
            std::visit([](auto& values) { values.shrink_to_fit(); }, m_values);
        }
        View view(std::size_t offset, std::size_t count) const {
            return std::visit(
                [&](const auto& values) { return View(values, offset, count); },
                m_values);
        }
        View view() const { return view(0, size()); }
        auto begin() const { return view().begin(); }
        auto end() const { return view().end(); }
        WideWeight operator[](std::size_t index) const { return view()[index]; }
        template <typename Function>
        decltype(auto) visit(Function&& function) const {
            return view().visit(std::forward<Function>(function));
        }
        void narrow(std::uint32_t row_nonzeros) {
            if (is_compact()) {
                return;
            }
            if (row_nonzeros <= 256) {
                narrow_to<ByteWeight>();
            } else if (row_nonzeros <= 65536) {
                narrow_to<ShortWeight>();
            }
        }

      private:
        template <typename NarrowWeight> void narrow_to() {
            const auto& wide = wide_values();
            std::vector<NarrowWeight> narrow;
            narrow.reserve(wide.size());
            for (const auto& value : wide) {
                narrow.emplace_back(value.row_inner_index(), value.weight());
            }
            // Assignment destroys the wide vector and its allocation only
            // after conversion succeeds. No construction-time slot is used.
            m_values = std::move(narrow);
        }
        std::variant<std::vector<WideWeight>, std::vector<ByteWeight>,
                     std::vector<ShortWeight>>
            m_values;
    };

} // namespace sasktran2::successive_orders
