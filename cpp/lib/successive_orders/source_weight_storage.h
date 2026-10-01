#pragma once

#include <algorithm>
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

    /** Finalized byte slots with exact IEEE weight representations.
     * The first record_count words hold records; the remaining words hold
     * escaped original double bits. No borrowed decoder pointer is retained.
     */
    struct EncodedSourceWeightStorage {
        std::vector<std::uint64_t> words;
        std::uint32_t record_count = 0;
        std::uint16_t base_exponent = 0;

        EncodedSourceWeightStorage() = default;
        EncodedSourceWeightStorage(const EncodedSourceWeightStorage&) = default;
        EncodedSourceWeightStorage&
        operator=(const EncodedSourceWeightStorage&) = default;
        EncodedSourceWeightStorage(EncodedSourceWeightStorage&& other) noexcept
            : words(std::move(other.words)),
              record_count(std::exchange(other.record_count, 0)),
              base_exponent(std::exchange(other.base_exponent, 0)) {}
        EncodedSourceWeightStorage&
        operator=(EncodedSourceWeightStorage&& other) noexcept {
            if (this != &other) {
                words = std::move(other.words);
                record_count = std::exchange(other.record_count, 0);
                base_exponent = std::exchange(other.base_exponent, 0);
            }
            return *this;
        }
        std::size_t size() const { return record_count; }
        std::size_t escape_count() const { return words.size() - record_count; }
        std::size_t capacity() const {
            return words.capacity() - escape_count();
        }
        void shrink_to_fit() { words.shrink_to_fit(); }
    };

    /** A decoder copied into each view, iterator and returned proxy. */
    struct SourceWeightDecoder {
        const std::uint64_t* escapes = nullptr;
        std::uint64_t base_shift = 0;

        std::uint64_t bits(std::uint64_t record) const {
            constexpr std::uint64_t payload_mask = (std::uint64_t{1} << 56) - 1;
            constexpr std::uint64_t escape_tag = std::uint64_t{15} << 52;
            const auto payload = record & payload_mask;
            if ((payload & escape_tag) == escape_tag) {
                return escapes[static_cast<std::uint32_t>(payload)];
            }
            return payload + base_shift;
        }
        bool operator==(const SourceWeightDecoder& other) const {
            return escapes == other.escapes && base_shift == other.base_shift;
        }
    };

    /** A value proxy; it never borrows an iterator or temporary view. */
    class EncodedSourceWeight {
      public:
        EncodedSourceWeight(std::uint64_t record, SourceWeightDecoder decoder)
            : m_record(record), m_decoder(decoder) {}
        std::uint32_t row_inner_index() const {
            return static_cast<std::uint32_t>(m_record >> 56);
        }
        std::uint64_t raw_weight_bits() const {
            return m_decoder.bits(m_record);
        }
        double weight() const {
            const auto bits = raw_weight_bits();
            double value;
            std::memcpy(&value, &bits, sizeof(value));
            return value;
        }

      private:
        std::uint64_t m_record;
        SourceWeightDecoder m_decoder;
    };

    /** A typed range for the encoded branch selected once per hot loop. */
    class EncodedSourceWeightArrayView {
      public:
        class Iterator {
          public:
            using iterator_category = std::input_iterator_tag;
            using value_type = EncodedSourceWeight;
            using difference_type = std::ptrdiff_t;
            using reference = EncodedSourceWeight;
            using pointer = void;

            EncodedSourceWeight operator*() const {
                return {*m_current, m_decoder};
            }
            Iterator& operator++() {
                ++m_current;
                return *this;
            }
            Iterator operator++(int) {
                auto previous = *this;
                ++*this;
                return previous;
            }
            bool operator==(const Iterator& other) const {
                return m_current == other.m_current &&
                       m_decoder == other.m_decoder;
            }
            bool operator!=(const Iterator& other) const {
                return !(*this == other);
            }

          private:
            friend class EncodedSourceWeightArrayView;
            Iterator(const std::uint64_t* current, SourceWeightDecoder decoder)
                : m_current(current), m_decoder(decoder) {}
            const std::uint64_t* m_current;
            SourceWeightDecoder m_decoder;
        };
        EncodedSourceWeightArrayView(const std::uint64_t* data,
                                     std::size_t size,
                                     SourceWeightDecoder decoder)
            : m_data(size == 0 ? nullptr : data), m_size(size),
              m_decoder(decoder) {}
        const std::uint64_t* data() const { return m_data; }
        std::size_t size() const { return m_size; }
        bool empty() const { return m_size == 0; }
        Iterator begin() const { return {m_data, m_decoder}; }
        Iterator end() const {
            return {m_size == 0 ? m_data : m_data + m_size, m_decoder};
        }
        EncodedSourceWeight operator[](std::size_t index) const {
            if (index >= m_size) {
                throw std::out_of_range("Encoded source weight view index");
            }
            return {m_data[index], m_decoder};
        }

      private:
        const std::uint64_t* m_data;
        std::size_t m_size;
        SourceWeightDecoder m_decoder;
    };

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
                if (m_encoded) {
                    const EncodedSourceWeight value(
                        static_cast<const std::uint64_t*>(m_data)[m_index],
                        m_decoder);
                    return {static_cast<int>(value.row_inner_index()),
                            value.weight()};
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
                       m_size == other.m_size &&
                       m_index_bytes == other.m_index_bytes &&
                       m_encoded == other.m_encoded &&
                       m_decoder == other.m_decoder;
            }
            bool operator!=(const Iterator& other) const {
                return !(*this == other);
            }

          private:
            friend class SourceWeightView;
            Iterator(const void* data, std::size_t size, std::size_t index,
                     unsigned char index_bytes, bool encoded,
                     SourceWeightDecoder decoder)
                : m_data(data), m_size(size), m_index(index),
                  m_index_bytes(index_bytes), m_encoded(encoded),
                  m_decoder(decoder) {}
            const void* m_data;
            std::size_t m_size;
            std::size_t m_index;
            unsigned char m_index_bytes;
            bool m_encoded;
            SourceWeightDecoder m_decoder;
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
        SourceWeightView(const EncodedSourceWeightStorage& values,
                         std::size_t offset, std::size_t size)
            : m_size(size), m_index_bytes(1), m_encoded(true) {
            if (values.record_count > values.words.size() ||
                values.base_exponent > 2032) {
                throw std::out_of_range("Encoded source weight record bounds");
            }
            if (offset > values.size() || size > values.size() - offset) {
                throw std::out_of_range("Encoded source weight view bounds");
            }
            m_data = size == 0 ? nullptr : values.words.data() + offset;
            // The escape table follows the WHOLE storage, never this subview.
            m_decoder = {values.words.empty()
                             ? nullptr
                             : values.words.data() + values.record_count,
                         static_cast<std::uint64_t>(values.base_exponent)
                             << 52};
        }
        std::size_t size() const { return m_size; }
        bool empty() const { return m_size == 0; }
        unsigned char index_bytes() const { return m_index_bytes; }
        bool is_encoded() const { return m_encoded; }
        Iterator begin() const {
            return {m_data, m_size, 0, m_index_bytes, m_encoded, m_decoder};
        }
        Iterator end() const {
            return {m_data,        m_size,    m_size,
                    m_index_bytes, m_encoded, m_decoder};
        }
        WideWeight operator[](std::size_t index) const {
            return *Iterator(m_data, m_size, index, m_index_bytes, m_encoded,
                             m_decoder);
        }
        template <typename Function>
        decltype(auto) visit(Function&& function) const {
            if (m_encoded) {
                return std::forward<Function>(function)(
                    EncodedSourceWeightArrayView(
                        static_cast<const std::uint64_t*>(m_data), m_size,
                        m_decoder));
            }
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
        bool m_encoded = false;
        SourceWeightDecoder m_decoder;
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
                    using Storage = std::decay_t<decltype(values)>;
                    if constexpr (std::is_same_v<Storage,
                                                 EncodedSourceWeightStorage>) {
                        return values.words.capacity() * sizeof(std::uint64_t);
                    } else {
                        return values.capacity() *
                               sizeof(typename Storage::value_type);
                    }
                },
                m_values);
        }
        std::size_t element_bytes() const {
            return std::visit(
                [](const auto& values) {
                    using Storage = std::decay_t<decltype(values)>;
                    if constexpr (std::is_same_v<Storage,
                                                 EncodedSourceWeightStorage>) {
                        return sizeof(std::uint64_t);
                    } else {
                        return sizeof(typename Storage::value_type);
                    }
                },
                m_values);
        }
        std::size_t payload_bytes() const {
            return size() * element_bytes() +
                   encoded_escape_count() * sizeof(std::uint64_t);
        }
        bool is_compact() const { return m_values.index() != 0; }
        static constexpr std::size_t encoded_header_growth_bytes() {
            using Legacy =
                std::variant<std::vector<WideWeight>, std::vector<ByteWeight>,
                             std::vector<ShortWeight>>;
            using Extended =
                std::variant<std::vector<WideWeight>, std::vector<ByteWeight>,
                             std::vector<ShortWeight>,
                             EncodedSourceWeightStorage>;
            return sizeof(Extended) - sizeof(Legacy);
        }
        bool is_encoded() const {
            return std::holds_alternative<EncodedSourceWeightStorage>(m_values);
        }
        std::size_t encoded_storage_bytes() const {
            const auto* encoded =
                std::get_if<EncodedSourceWeightStorage>(&m_values);
            return encoded == nullptr
                       ? 0
                       : encoded->words.capacity() * sizeof(std::uint64_t);
        }
        std::size_t encoded_escape_count() const {
            const auto* encoded =
                std::get_if<EncodedSourceWeightStorage>(&m_values);
            return encoded == nullptr ? 0 : encoded->escape_count();
        }
        std::uint16_t encoded_base_exponent() const {
            const auto* encoded =
                std::get_if<EncodedSourceWeightStorage>(&m_values);
            return encoded == nullptr ? 0 : encoded->base_exponent;
        }
        std::uint32_t encoded_record_count() const {
            const auto* encoded =
                std::get_if<EncodedSourceWeightStorage>(&m_values);
            return encoded == nullptr ? 0 : encoded->record_count;
        }
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
                if (!try_encode()) {
                    narrow_to<ByteWeight>();
                }
            } else if (row_nonzeros <= 65536) {
                narrow_to<ShortWeight>();
            }
        }

      private:
        static std::uint64_t raw_bits(const WideWeight& value) {
            const double weight = value.weight();
            std::uint64_t bits;
            std::memcpy(&bits, &weight, sizeof(bits));
            return bits;
        }
        static bool window_eligible(std::uint64_t bits) {
            // Positive finite nonzero values only; this is integer
            // classification. Signed zeros, negative values and specials
            // escape with every original bit unchanged.
            return bits != 0 && bits < 0x7ff0000000000000ULL;
        }
        bool try_encode() {
            if constexpr (!std::numeric_limits<double>::is_iec559 ||
                          std::numeric_limits<double>::digits != 53 ||
                          std::numeric_limits<double>::max_exponent != 1024) {
                return false;
            }
            const auto& wide = wide_values();
            const auto count = wide.size();
            // Short arrays cannot repay the larger container header.
            // The byte-slot fallback still validates every slot.
            if (count <= encoded_header_growth_bytes() ||
                count > std::numeric_limits<std::uint32_t>::max()) {
                return false;
            }
            // Only touched histogram bins are initialized. Most individual
            // rays occupy a small exponent range; scanning all 2047 bins for
            // each short ray would add needless construction work.
            std::array<std::uint32_t, 2047> histogram;
            std::array<std::uint64_t, 32> present{};
            unsigned int minimum = 2047;
            unsigned int maximum = 0;
            for (const auto& value : wide) {
                if (value.row_inner_index() > 255) {
                    throw std::out_of_range(
                        "Source interpolation slot overflow");
                }
                const auto bits = raw_bits(value);
                if (!window_eligible(bits)) {
                    continue;
                }
                const auto exponent = static_cast<unsigned int>(bits >> 52);
                const auto mask = std::uint64_t{1} << (exponent & 63);
                auto& visited = present[exponent >> 6];
                if ((visited & mask) == 0) {
                    visited |= mask;
                    histogram[exponent] = 0;
                }
                ++histogram[exponent];
                minimum = std::min(minimum, exponent);
                maximum = std::max(maximum, exponent);
            }
            if (minimum == 2047) {
                return false;
            }
            auto bin = [&](unsigned int exponent) -> std::uint64_t {
                const auto mask = std::uint64_t{1} << (exponent & 63);
                return (present[exponent >> 6] & mask) == 0
                           ? 0
                           : histogram[exponent];
            };
            const auto first_base = minimum > 14 ? minimum - 14 : 0;
            const auto last_base = std::min(maximum, 2032U);
            std::uint64_t window_count = 0;
            for (unsigned int exponent = first_base; exponent < first_base + 15;
                 ++exponent) {
                window_count += bin(exponent);
            }
            auto best_base = first_base;
            auto best_count = window_count;
            for (unsigned int base = first_base + 1; base <= last_base;
                 ++base) {
                window_count -= bin(base - 1);
                window_count += bin(base + 14);
                if (window_count > best_count) {
                    best_count = window_count;
                    best_base = base;
                }
            }
            const auto escapes = count - static_cast<std::size_t>(best_count);
            // Exact new payload plus the 8-byte container-header delta must
            // beat ordinary byte-slot records before allocating anything.
            const auto old_payload = static_cast<std::uint64_t>(count) * 9;
            const auto new_payload =
                static_cast<std::uint64_t>(count + escapes) * 8;
            if (new_payload + encoded_header_growth_bytes() >= old_payload) {
                return false;
            }
            EncodedSourceWeightStorage encoded;
            encoded.record_count = static_cast<std::uint32_t>(count);
            encoded.base_exponent = static_cast<std::uint16_t>(best_base);
            encoded.words.resize(count + escapes);
            constexpr std::uint64_t mantissa_mask =
                (std::uint64_t{1} << 52) - 1;
            std::uint32_t escape_index = 0;
            for (std::size_t index = 0; index < count; ++index) {
                const auto& value = wide[index];
                const auto bits = raw_bits(value);
                const auto exponent =
                    static_cast<unsigned int>((bits >> 52) & 0x7ff);
                std::uint64_t payload;
                if (window_eligible(bits) && exponent >= best_base &&
                    exponent < best_base + 15) {
                    payload = (bits & mantissa_mask) |
                              (static_cast<std::uint64_t>(exponent - best_base)
                               << 52);
                } else {
                    payload = (std::uint64_t{15} << 52) | escape_index;
                    encoded.words[count + escape_index] = bits;
                    ++escape_index;
                }
                encoded.words[index] =
                    payload |
                    (static_cast<std::uint64_t>(value.row_inner_index()) << 56);
            }
            if (escape_index != escapes) {
                throw std::logic_error("Source weight codec escape mismatch");
            }
            // Honor real allocator capacity, not just requested element count.
            const auto allocated =
                static_cast<std::uint64_t>(encoded.words.capacity()) * 8;
            if (allocated + encoded_header_growth_bytes() >= old_payload) {
                return false;
            }
            // This is the only mutation of the existing variant. Failed
            // allocation/conversion leaves the original wide records intact.
            m_values = std::move(encoded);
            return true;
        }
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
                     std::vector<ShortWeight>, EncodedSourceWeightStorage>
            m_values;
    };

} // namespace sasktran2::successive_orders
